#!/usr/bin/env python3
"""Convert lya2pcf data*.npy dictionaries to cTreeBalls lya-ascii."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys
from typing import Iterable

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Flatten lya2pcf forest objects into x y z delta weight forest_id "
            "rows for the cTreeBalls Ly-alpha addon."
        )
    )
    parser.add_argument("inputs", nargs="+", type=Path, help="lya2pcf data*.npy files")
    parser.add_argument("-o", "--output", required=True, type=Path)
    parser.add_argument(
        "--lya2pcf-source",
        type=Path,
        help="directory containing forest_class.py, needed while unpickling",
    )
    return parser.parse_args()


def forest_key(forest: object) -> tuple[str, str]:
    name = getattr(forest, "name", None)
    if name is None:
        raise ValueError("a forest object has no 'name' attribute")
    return type(name).__name__, repr(name)


def forests_in_file(path: Path) -> Iterable[object]:
    loaded = np.load(path, allow_pickle=True)
    try:
        data = loaded.item()
    except ValueError as exc:
        raise ValueError(f"{path}: expected a saved dictionary") from exc
    if not isinstance(data, dict):
        raise ValueError(f"{path}: expected a dictionary, got {type(data).__name__}")
    for forests in data.values():
        for forest in forests:
            yield forest


def forest_rows(forest: object, forest_id: int) -> Iterable[tuple[float, ...]]:
    required = ("rx", "ry", "rz", "we", "dw")
    arrays = {}
    for name in required:
        if not hasattr(forest, name):
            raise ValueError(f"forest {getattr(forest, 'name', '?')!r} lacks {name!r}")
        arrays[name] = np.asarray(getattr(forest, name), dtype=np.float64).reshape(-1)
    lengths = {array.size for array in arrays.values()}
    if len(lengths) != 1:
        raise ValueError(f"forest {getattr(forest, 'name', '?')!r} has unequal array lengths")

    for index, (x, y, z, weight, weighted_delta) in enumerate(
        zip(arrays["rx"], arrays["ry"], arrays["rz"], arrays["we"], arrays["dw"])
    ):
        values = np.array((x, y, z, weight, weighted_delta), dtype=np.float64)
        if not np.all(np.isfinite(values)):
            raise ValueError(
                f"forest {getattr(forest, 'name', '?')!r}, pixel {index}: non-finite value"
            )
        if weight < 0.0:
            raise ValueError(
                f"forest {getattr(forest, 'name', '?')!r}, pixel {index}: negative weight"
            )
        if np.linalg.norm((x, y, z)) <= 0.0:
            raise ValueError(
                f"forest {getattr(forest, 'name', '?')!r}, pixel {index}: zero distance"
            )
        if weight == 0.0:
            if weighted_delta != 0.0:
                raise ValueError(
                    f"forest {getattr(forest, 'name', '?')!r}, pixel {index}: "
                    "zero weight with nonzero dw"
                )
            delta = 0.0
        else:
            delta = weighted_delta / weight
        yield x, y, z, delta, weight, forest_id


def main() -> int:
    args = parse_args()
    if args.lya2pcf_source is not None:
        source = args.lya2pcf_source.expanduser().resolve()
        if not (source / "forest_class.py").is_file():
            raise SystemExit(f"{source} does not contain forest_class.py")
        sys.path.insert(0, os.fspath(source))

    output = args.output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    forest_ids: dict[tuple[str, str], int] = {}
    row_count = 0

    try:
        with temporary.open("w", encoding="ascii") as stream:
            stream.write("# x y z delta weight forest_id\n")
            for input_path in args.inputs:
                for forest in forests_in_file(input_path.expanduser().resolve()):
                    key = forest_key(forest)
                    forest_id = forest_ids.setdefault(key, len(forest_ids))
                    for row in forest_rows(forest, forest_id):
                        stream.write(
                            "{:.17g} {:.17g} {:.17g} {:.17g} {:.17g} {}\n".format(*row)
                        )
                        row_count += 1
        if row_count == 0:
            raise ValueError("the input dictionaries contain no forest pixels")
        temporary.replace(output)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise

    print(f"wrote {row_count} pixels from {len(forest_ids)} forests to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Interactive terminal interface for the offline and optional local-AI tutor."""
import argparse
import json
import math
import os
from pathlib import Path
import re
import shlex
import sys
import sqlite3
import struct
import subprocess
import zipfile

from .catalog import inspect_catalog, read_run
from .knowledge import Knowledge
from .local_model import LocalModel
from .session import Answer, Tutor

HELP = """Ask a question in ordinary English, or use:
  /topics                  List the reviewed scientific explanations
  /topic NAME              Read one explanation
  /search WORDS            Search current documentation and README files
  /sources                 List sources for the last answer
  /source NUMBER           Show a source excerpt from the last answer
  /context key=value ...   Set field, geometry, statistic, catalog, points,
                           threads or machine (quoted values are accepted)
  /context                 Show the current scientific context
  /recommend               Compare compatible engine candidates
  /engines [shear|scalar|lya]  List methods and check the compiled executable
  /defaults                Explain the approximation defaults
  /build                   Show relevant settings from this cballs executable
  /catalog PATH            Inspect FITS/NPZ headers, without loading the arrays
  /run PATH                Explain a saved approximation-penalty JSON or summary
  /reindex                 Refresh the local documentation database
  /save PATH               Save this session as JSON (refuses to overwrite)
  /reset                   Clear the conversation and scientific context
  /help                    Show these commands
  /quit                    Exit; Ctrl-D also exits

Example: /context field=shear geometry=sphere statistic=3pcf threads=16
Offline mode retrieves documentation and reviewed explanations; it is not a
generative language model. It never executes suggested analysis commands.
"""


def clean_terminal(text):
    return re.sub(r"[\x00-\x08\x0b-\x1f\x7f]", "", str(text))


def render(answer):
    output = answer.text
    if answer.notice:
        output = answer.notice+"\n\n"+output
    if answer.sources:
        output += "\n\nSources:\n"+"\n".join(f"  [{i}] {source}" for i, source in enumerate(answer.sources, 1))
    return clean_terminal(output)


def save_session(tutor, filename):
    path = Path(filename).expanduser()
    payload = dict(schema_version=1, context=tutor.context, turns=tutor.history,
                   source_root=str(tutor.root), mode="local-model" if tutor.model else "offline")
    with path.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    return f"Saved session to {path.resolve()}"


def run_report(tutor, filename):
    path, data = read_run(filename)
    lines = []
    if "observables" in data and "target_percent" in data:
        if (not isinstance(data["observables"], dict) or
                not isinstance(data["target_percent"], (int, float)) or
                not math.isfinite(data["target_percent"]) or data["target_percent"] < 0):
            raise ValueError("Invalid approximation-penalty report: expected observables and a finite target.")
        lines.append(f"Approximation penalty: {data.get('status', 'unknown')}; target {data['target_percent']}%")
        if not data["observables"]:
            lines.append("No measured reference comparison is stored here. UNMEASURED does not mean zero error.")
        for name, row in data["observables"].items():
            if not isinstance(row, dict):
                raise ValueError("Invalid observable in approximation-penalty report.")
            percent = row.get("relative_l2_percent")
            valid = isinstance(percent, (int, float)) and math.isfinite(percent) and percent >= 0
            value = f"{percent:.6g}%" if valid else "undefined"
            state = "within target" if valid and percent <= data["target_percent"] and row.get("relative_target_met") else "above target or undefined"
            lines.append(f"  {name}: relative L2 {value}; {state}; absolute L2 {row.get('absolute_l2')}")
        if "native_wall_speedup" in data:
            lines.append(f"Recorded native wall speedup: {data['native_wall_speedup']:.4g}x")
        lines += ["", "This evidence applies only to the matched catalogue, estimator, binning and retained multipoles. It does not establish a per-bin or full-catalogue error bound."]
    else:
        metadata = data.get("metadata", data)
        if not isinstance(metadata, dict):
            raise ValueError("Run metadata must be a JSON object.")
        if "engine" in metadata:
            lines.append("Engine: "+str(metadata["engine"]))
        for key in ("opening_tolerance", "effective_smoothing", "shear_hierarchical_reuse",
                    "scalar_hierarchical_reuse", "timings", "qualification"):
            if key in metadata:
                item = metadata[key]
                if key == "qualification" and isinstance(item, dict):
                    item = {k: item[k] for k in ("status", "tolerances", "limitation") if k in item}
                lines.append(key+": "+json.dumps(item, ensure_ascii=False)[:2500])
        if "engines" in data:
            lines.append("Engines in summary: "+", ".join(data["engines"]))
        if not lines:
            raise ValueError("This JSON has no recognized run information or approximation report.")
        lines.append("Use /run on ENGINE/qualification/approximation-penalty.json for the per-observable error table.")
    return Answer("Saved run: "+str(path), "\n".join(lines), [str(path)], context=dict(tutor.context))


def command(tutor, text):
    name, _, argument = text.strip().partition(" ")
    argument = argument.strip()
    if name in {"/quit", "/exit"}:
        return None
    if name == "/help":
        return HELP
    if name == "/topics":
        return "\n".join(f"  {p.topic}: {p.title}" for p in tutor.knowledge.passages if p.kind == "answer")
    if name == "/topic":
        answer = tutor.topic(argument)
    elif name == "/defaults":
        answer = tutor.topic("defaults")
    elif name == "/recommend":
        answer = tutor.recommend()
    elif name == "/run":
        answer = run_report(tutor, path_argument(argument))
    elif name == "/sources":
        return "\n".join(f"[{i}] {s}" for i, s in enumerate(tutor.last.sources, 1)) if tutor.last else "No answer yet."
    elif name == "/source":
        if not tutor.last or not argument.isdigit() or not 1 <= int(argument) <= len(tutor.last.sources):
            raise ValueError("Choose a source number from /sources.")
        source = tutor.last.sources[int(argument)-1]
        filename, _, number = source.rpartition(":")
        if not number.isdigit():
            filename, number = source, "1"
        path = (tutor.root/filename).resolve()
        if not path.is_relative_to(tutor.root) or not path.is_file() or path.stat().st_size > 2_000_000:
            raise ValueError("This source is outside the indexed checkout; open its cited file directly.")
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        start = max(int(number)-1, 0)
        return source+"\n"+"\n".join(f"{i+1}: {line}" for i, line in enumerate(lines[start:start+70], start))
    elif name == "/search":
        if not argument:
            raise ValueError("Use /search followed by an engine, parameter or scientific term.")
        matches = tutor.knowledge.search(argument)
        answer = Answer(argument, "\n\n".join(f"{p.title}\n{p.text[:1000]}" for _, _, p in matches)
                        or "No matching local passages.", [p.citation for _, _, p in matches])
    elif name == "/context":
        if argument:
            tutor.set_context(argument)
        return json.dumps(tutor.context, indent=2, ensure_ascii=False) if tutor.context else "No context set. Try /context field=shear geometry=sphere statistic=3pcf"
    elif name == "/engines":
        return tutor.engines(argument)
    elif name == "/build":
        output = tutor.probe("make-info")
        return "\n".join(line for line in output.splitlines() if
            re.match(r"\s*(THETA|CBALLS_(SHEAR|SCALAR)_|OPENMP|MPI|SINGLEPON|SMOOTHPIVOTON|BALLS4SCANLEVON)\w*\s*=", line)
            or any(phrase in line for phrase in ("runtime theta defaults", "reuse compiled defaults"))) or output[:8000]
    elif name == "/catalog":
        info = inspect_catalog(path_argument(argument))
        tutor.context["catalog"] = info["path"]
        for key in ("points", "field", "geometry"):
            tutor.context.pop(key, None)
        if info["rows"]:
            tutor.context["points"] = info["rows"]
        if len(info["fields"]) == 1:
            tutor.context["field"] = info["fields"][0]
        if info["geometry_hint"]:
            tutor.context["geometry"] = info["geometry_hint"]
        return json.dumps(info, indent=2, ensure_ascii=False)+"\nSelect the intended field and geometry with /context, then use /recommend."
    elif name == "/reset":
        tutor.context.clear(); tutor.history.clear(); tutor.last = None
        return "Conversation and scientific context cleared."
    elif name == "/reindex":
        root, database = tutor.root, tutor.knowledge.database
        tutor.knowledge.close()
        tutor.knowledge = Knowledge(root, database, rebuild=True)
        tutor._available = None
        return f"Indexed {len(tutor.knowledge.file_hashes)} source files and {len(tutor.knowledge.passages)} passages."
    elif name == "/save":
        return save_session(tutor, path_argument(argument))
    else:
        raise ValueError("Unknown command. Use /help, or ask a question without the leading slash.")
    tutor.last = answer
    tutor.history.append(dict(question=text, answer=answer.text, sources=answer.sources, mode=answer.mode))
    return render(answer)


def path_argument(argument):
    items = shlex.split(argument)
    if len(items) != 1:
        raise ValueError('Provide one path; enclose a path containing spaces in quotes.')
    return items[0]


def parser():
    p = argparse.ArgumentParser(description="cTreeBalls IA tutor: offline documentation and reviewed scientific explanations, with optional local Ollama synthesis.")
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--cballs", type=Path, help="native executable to inspect for available engines")
    p.add_argument("--ask", help="answer one question and exit")
    p.add_argument("--json", action="store_true", help="machine-readable one-question answer")
    p.add_argument("--database", default=os.environ.get("CBALLS_TUTOR_DATABASE"), help="SQLite cache path or :memory:; default ROOT/build/ia-tutor/knowledge.sqlite3")
    p.add_argument("--reindex", action="store_true", help="rebuild the database and exit unless --ask is supplied")
    p.add_argument("--backend", choices=("offline", "ollama"), default=os.environ.get("CBALLS_TUTOR_BACKEND", "offline"))
    p.add_argument("--model", default=os.environ.get("CBALLS_TUTOR_MODEL"), help="name of an already installed local model")
    p.add_argument("--model-url", default=os.environ.get("CBALLS_TUTOR_OLLAMA_URL", "http://127.0.0.1:11434"))
    p.add_argument("--timeout", type=float, default=60, help="local-model request timeout in seconds")
    p.add_argument("--context", action="append", default=[], help='initial context, e.g. "field=shear geometry=sphere statistic=3pcf"')
    p.add_argument("--transcript", type=Path, help="explicitly save session JSON on exit; refuses to overwrite")
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    if args.json and not args.ask:
        parser().error("--json requires --ask")
    if not (args.root/"tutor/knowledge/questions.json").is_file():
        parser().error("--root must point to a source checkout containing tutor/knowledge/questions.json")
    for key in ("OMPI_COMM_WORLD_SIZE", "PMI_SIZE", "MV2_COMM_WORLD_SIZE"):
        if os.environ.get(key, "1").isdigit() and int(os.environ.get(key, "1")) > 1:
            parser().error("The tutor is a serial conversation; launch it outside mpiexec/srun.")
    if (os.environ.get("SLURM_PROCID") is not None and os.environ.get("SLURM_NTASKS", "1").isdigit()
            and int(os.environ.get("SLURM_NTASKS", "1")) > 1):
        parser().error("The tutor needs one terminal, not multiple Slurm tasks.")
    tutor = None
    try:
        if args.transcript and args.transcript.exists():
            raise ValueError("Transcript already exists; choose a new filename.")
        model = LocalModel(args.model, args.model_url, args.timeout) if args.backend == "ollama" else None
        knowledge = Knowledge(args.root, args.database, args.reindex)
        tutor = Tutor(knowledge, args.cballs, model)
        for item in args.context:
            tutor.set_context(item)
        if args.ask:
            answer = tutor.ask(args.ask)
            print(json.dumps(answer.as_dict(), indent=2, ensure_ascii=False, allow_nan=False) if args.json else render(answer))
        elif args.reindex:
            print(f"Indexed {len(knowledge.file_hashes)} source files and {len(knowledge.passages)} passages in {knowledge.database}")
        else:
            print("cTreeBalls IA tutor — "+("offline documentation mode" if model is None else "local Ollama model: "+args.model))
            print("Ask a scientific question. /help lists commands; /quit exits.")
            print(f"Knowledge: {len(knowledge.file_hashes)} source files; {len(knowledge.passages)} passages. Questions are not saved unless you use /save or --transcript.")
            if knowledge.cache_notice:
                print(knowledge.cache_notice)
            while True:
                if sys.stdin.isatty():
                    try:
                        text = input("You> ")
                    except EOFError:
                        break
                else:
                    text = sys.stdin.readline()
                    if not text:
                        break
                text = text.strip()
                if not text:
                    continue
                if text.lower() in {"quit", "exit", "/quit", "/exit"}:
                    break
                try:
                    result = command(tutor, text) if text.startswith("/") else render(tutor.ask(text))
                    if result is None:
                        break
                    print("\n"+clean_terminal(result)+"\n")
                except (ValueError, OSError, KeyError, TypeError, sqlite3.Error, struct.error, zipfile.BadZipFile, subprocess.SubprocessError) as exc:
                    print("Tutor: "+clean_terminal(exc), file=sys.stderr)
        if args.transcript:
            message = save_session(tutor, args.transcript)
            if not args.json:
                print(message)
        return 0
    except KeyboardInterrupt:
        print("\nTutor interrupted. No analysis job was started.", file=sys.stderr)
        return 130
    except (ValueError, OSError, KeyError, TypeError, sqlite3.Error, struct.error, zipfile.BadZipFile, subprocess.SubprocessError) as exc:
        print("Tutor: "+clean_terminal(exc), file=sys.stderr)
        return 2
    finally:
        if tutor:
            tutor.knowledge.close()


if __name__ == "__main__":
    raise SystemExit(main())

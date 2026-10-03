"""Shared elapsed-time and process-CPU accounting for all native drivers.

Wall time is a critical-path measurement, CPU is consumed process time (including
OpenMP workers). MPI CPU totals sum participating ranks. Catalog loading,
registration, plotting and launcher overhead are outside the timed Run scope.
"""
import math


def _seconds(value, name):
    value = float(value)
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"{name} must be finite nonnegative seconds")
    return value


def timing_metadata(setup_wall, setup_cpu, compute_wall, compute_cpu, scope):
    values = dict(zip(("setup_wall_time", "setup_cpu_time", "compute_wall_time",
                       "compute_cpu_time"), (setup_wall, setup_cpu, compute_wall, compute_cpu)))
    result = {name: _seconds(value, name) for name, value in values.items()}
    result.update(total_wall_time=result["setup_wall_time"]+result["compute_wall_time"],
                  total_cpu_time=result["setup_cpu_time"]+result["compute_cpu_time"],
                  timing_scope=scope)
    return result


def aggregate_rank_timings(rows):
    """Keep the maximum per-rank total, not the sum of separate maxima."""
    if not rows:
        raise ValueError("at least one participating rank is required")
    ranks = [row.get("rank", index) for index, row in enumerate(rows)]
    if len(set(ranks)) != len(rows):
        raise ValueError("each participating rank must occur exactly once")
    required = ("setup_wall_time", "setup_cpu_time", "compute_wall_time",
                "compute_cpu_time", "total_wall_time", "total_cpu_time",
                "native_reported_cpu_time")
    result = {}
    for name in (*required, "native_mainloop_wall_time", "native_mainloop_cpu_time"):
        if name not in required and not any(name in row for row in rows):
            continue
        if not all(name in row for row in rows):
            raise ValueError(f"missing {name} on a participating rank")
        values = [_seconds(row[name], name) for row in rows]
        result[name] = (max if "wall" in name else sum)(values)
    scopes = {row.get("timing_scope", "native setup and Python Run") for row in rows}
    if len(scopes) != 1:
        raise ValueError("participating ranks must use the same timing scope")
    result.update(ranks=len(rows), rank_timings=list(rows),
                  timing_scope=scopes.pop()+
                  "; wall=max(participating ranks), CPU=sum(participating ranks)")
    return result

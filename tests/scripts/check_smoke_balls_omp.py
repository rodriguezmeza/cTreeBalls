# tests/scripts/check_smoke_balls_omp.py
import sys
from pathlib import Path
import numpy as np

out = Path(sys.argv[1])

required = [
    out / "histNN.txt",
    out / "histXi2pcf.txt",
]

for path in required:
    if not path.exists():
        raise SystemExit(f"missing output file: {path}")

hist_nn = np.loadtxt(out / "histNN.txt")
hist_xi = np.loadtxt(out / "histXi2pcf.txt")

if hist_nn.size == 0:
    raise SystemExit("histNN.txt is empty")

if not np.isfinite(hist_nn).all():
    raise SystemExit("histNN.txt contains non-finite values")

# Usually count columns after radius/bin column. Adjust if your format differs.
count_sum = np.sum(hist_nn[:, 1:]) if hist_nn.ndim == 2 and hist_nn.shape[1] > 1 else np.sum(hist_nn)

if count_sum <= 0:
    raise SystemExit(f"histNN has no positive counts: sum={count_sum}")

if not np.isfinite(hist_xi).all():
    raise SystemExit("histXi2pcf.txt contains non-finite values")

print(f"balls-omp smoke test passed: count_sum={count_sum:g}")

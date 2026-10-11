#!/usr/bin/env python3
"""Run the source-checkout cTreeBalls tutor without importing cyballs."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tutor.cli import main

if __name__ == "__main__":
    raise SystemExit(main())

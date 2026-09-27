#!/usr/bin/env python3
"""Print a checked per-rank forecast without allocating a catalog."""
import argparse,json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from cyballs import resource_plan
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--engine',required=True);p.add_argument('--pixels',required=True,type=int)
p.add_argument('--threads',type=int,default=1);p.add_argument('--parameters',type=Path)
a=p.parse_args();params=json.loads(a.parameters.read_text()) if a.parameters else {}
print(json.dumps(resource_plan(a.engine,a.pixels,a.threads,params),indent=2))

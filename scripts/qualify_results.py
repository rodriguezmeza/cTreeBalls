#!/usr/bin/env python3
"""Qualify saved user-driver results against a matching exact reference run.

The result packet contains the actual arrays, geometry and input fingerprints.
Native CLI metadata alone is deliberately insufficient evidence.
"""
import argparse
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from cyballs import load_result_packet,qualification_report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidate',required=True,type=Path)
    p.add_argument('--reference',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--rtol',type=float,default=.02);p.add_argument('--atol',type=float,default=1e-10)
    a=p.parse_args()
    report=qualification_report(load_result_packet(a.candidate),load_result_packet(a.reference),rtol=a.rtol,atol=a.atol)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(report['status'])
    for name,row in report['observables'].items():
        print(f"{name}: {row['status']}; relative L2={row['relative_l2']}; weak bins={row['weak_reference_bins']}; pointwise exceedances={row['bins_exceeding_pointwise_tolerance']}")
    return 0 if report['status']=='QUALIFIED' else 1
if __name__=='__main__':raise SystemExit(main())

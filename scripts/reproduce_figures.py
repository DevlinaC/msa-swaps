#!/usr/bin/env python3
"""Regenerate available manuscript plots from frozen local tables and renders."""
import argparse
import os
from pathlib import Path
import subprocess
import sys
ROOT=Path(__file__).resolve().parents[1]
COMMANDS={
    'S1':[['scripts/plot_s1.py']],
    'S2':[['scripts/plot_s2_s3.py','--figure','2']],
    'S3':[['scripts/plot_s2_s3.py','--figure','3']],
    'S4':[['scripts/plot_s4.py']],
    'S6':[['workflows/pairformer/scripts/plot_summary.py','--root','workflows/pairformer']],
    'S9':[['workflows/chimera/scripts/export_panels.py','--panel','all'],
          ['workflows/chimera/scripts/plot_chimera_s7.py']],
}
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--figure',choices=['all']+list(COMMANDS),default='all')
    args=parser.parse_args()
    env=dict(os.environ)
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
        env[key]='1'
    for name in COMMANDS if args.figure=='all' else [args.figure]:
        for command in COMMANDS[name]:
            print('Generating',name,':', ' '.join(command),flush=True)
            subprocess.run([sys.executable]+command,cwd=ROOT,env=env,check=True)
if __name__=='__main__':main()

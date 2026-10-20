#!/usr/bin/env python3
"""Create a Markdown performance report from soil-shell benchmark CSV data."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path
import statistics

LABELS={
    'analytical':'Analytical / semi-analytical approximation',
    'marching':'Steady spatial marching',
    'newton':'Coupled steady Newton',
    'theta_relaxation':'Implicit θ pseudo-time relaxation',
}
REQUIRED={'method','repeat','seconds','steps','newton_iterations',
          'scaled_steady_residual','p_out_Pa','Tg_out_K'}

def read_runs(path):
    groups=defaultdict(list)
    with path.open(newline='',encoding='utf-8') as stream:
        reader=csv.DictReader(stream)
        if not REQUIRED.issubset(reader.fieldnames or ()):
            raise ValueError(f'{path} is not soil_shell_benchmark CSV data')
        for row in reader:
            method=row['method']
            if method not in LABELS:
                raise ValueError(f'unknown soil-shell method: {method}')
            groups[method].append({
                'repeat':int(row['repeat']),'seconds':float(row['seconds']),
                'steps':int(row['steps']),'iterations':int(row['newton_iterations']),
                'residual':float(row['scaled_steady_residual']),
                'pressure':float(row['p_out_Pa']),'temperature':float(row['Tg_out_K'])})
    if not groups:
        raise ValueError(f'no timing rows found in {path}')
    return groups

def duration(seconds):
    if seconds>=1: return f'{seconds:.4f} s'
    if seconds>=1e-3: return f'{seconds*1e3:.3f} ms'
    if seconds>=1e-6: return f'{seconds*1e6:.3f} μs'
    return f'{seconds*1e9:.3f} ns'

def make_report(groups,csv_path):
    medians={k:statistics.median(r['seconds'] for r in v) for k,v in groups.items()}
    fastest=min(medians.values())
    lines=[
        '# Fully thermal soil-shell initializer performance report','',
        f'Source data: `{csv_path.name}`. Values summarize repeated runs; timing is measured on the machine that ran the benchmark.','',
        'Default case: 100 km pipe, 16 axial cells, 0.6 m diameter, 8 MPa inlet pressure, 288.15 K inlet gas, 278.15 K far-field soil, 100 kg/s outlet flow, four soil shells, flat elevation. Newton uses the analytical profile as its initial guess. θ relaxation uses backward Euler (θ=1) and the configured pseudo-time step.','',
        '| Method | Runs | Median | Mean | Std. dev. | Min–max | Work | Max residual | Outlet p / Tg | vs. fastest |',
        '|---|---:|---:|---:|---:|---:|---|---:|---|---:|']
    for method,label in LABELS.items():
        if method not in groups: continue
        runs=groups[method]; times=[r['seconds'] for r in runs]
        work=sorted({(r['iterations'],r['steps']) for r in runs})
        if method=='theta_relaxation':
            work_text='; '.join(f'{steps} steps, {iterations} inner Newton iterations' for iterations,steps in work)
        elif method=='newton':
            work_text='; '.join(f'{iterations} Newton iterations' for iterations,_ in work)
        else: work_text='direct evaluation'
        p=statistics.median(r['pressure'] for r in runs)/1e6
        t=statistics.median(r['temperature'] for r in runs)
        sd=statistics.stdev(times) if len(times)>1 else 0.0
        lines.append(f'| {label} | {len(runs)} | {duration(statistics.median(times))} | {duration(statistics.mean(times))} | {duration(sd)} | {duration(min(times))}–{duration(max(times))} | {work_text} | {max(r["residual"] for r in runs):.3e} | {p:.6f} MPa / {t:.4f} K | {medians[method]/fastest:.2f}× |')
    lines.extend(['','Timing ratios compare medians against the fastest method. The analytical approximation is a continuum estimate and is not expected to satisfy the discrete steady residual as tightly as the discrete solvers. Runtime is environment-dependent; do not compare figures across different hardware or compiler settings as if they were controlled benchmarks.',''])
    return '\n'.join(lines)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('csv',type=Path,help='CSV emitted by soil_shell_benchmark')
    parser.add_argument('-o','--output',type=Path,help='Markdown path; defaults to <csv stem>_report.md')
    args=parser.parse_args(); output=args.output or args.csv.with_name(f'{args.csv.stem}_report.md')
    report=make_report(read_runs(args.csv),args.csv)
    output.write_text(report,encoding='utf-8')
    print(report); print(f'\nWrote {output}')

if __name__=='__main__': main()

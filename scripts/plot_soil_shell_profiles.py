#!/usr/bin/env python3
# Copyright (C) 2026 Richard Myers
# SPDX-License-Identifier: AGPL-3.0-only
#
# Licensed under the GNU Affero General Public License, version 3 only.
# See the LICENSE file for the complete terms.
#
# A separate commercial license for proprietary use may be available by
# written agreement with the relevant rights holders. This notice does not
# grant that license, guarantee its availability, or set fees, royalties,
# or revenue shares.

"""Plot pressure, mass flow, gas/wall/soil temperatures and pipe elevation."""
import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MaxNLocator

METHODS=(('Analytical / semi-analytical approximation','analytical'),
         ('Steady spatial marching','marching'),
         ('Coupled steady Newton','newton'),
         ('Implicit θ pseudo-time relaxation','relaxation'))
COLORS={'p':'#1769aa','m':'#c44e52','tg':'#dd8b00','tw':'#7b61a8','z':'#477553'}
SOIL_COLORS=('#168c88','#6a9f45','#a77b36','#ab5b76','#6877aa','#8a6b52')

def read_profiles(path):
    with path.open(newline='',encoding='utf-8') as stream:
        rows=list(csv.DictReader(stream))
    if not rows or 'x_node_m' not in rows[0]:
        raise ValueError(f'no profile data found in {path}')
    return rows

def plot_profiles(rows,output):
    nodes=[r for r in rows if r['x_node_m']]
    faces=[r for r in rows if r['x_face_m']]
    x=[float(r['x_node_m'])/1000 for r in nodes]
    xf=[float(r['x_face_m'])/1000 for r in faces]
    soil_columns=[k for k in rows[0] if k.startswith('Ts0_') and k.endswith('_K')]
    n_shells=len(soil_columns)
    fig,axes=plt.subplots(5,1,figsize=(13.5,17.0),sharex=True,
        gridspec_kw={'height_ratios':[1,1,1,1,.62]})
    flows=[float(r[f'm_{key}_kg_s']) for _,key in METHODS for r in faces]
    center=sum(flows)/len(flows); pad=max(1.25*max(abs(v-center) for v in flows),abs(center)*1e-4,1e-6)
    legend_handles=None
    for ax,(title,key) in zip(axes[:4],METHODS):
        maxis=ax.twinx(); taxis=ax.twinx(); taxis.spines['right'].set_position(('outward',94))
        p,=ax.plot(x,[float(r[f'p_{key}_Pa'])/1e6 for r in nodes],color=COLORS['p'],lw=1.8,label='Pressure')
        m,=maxis.plot(xf,[float(r[f'm_{key}_kg_s']) for r in faces],color=COLORS['m'],ls='--',marker='o',ms=2,lw=1.4,label='Mass flow')
        tg,=taxis.plot(x,[float(r[f'Tg_{key}_K']) for r in nodes],color=COLORS['tg'],lw=1.7,label='Gas temperature')
        tw,=taxis.plot(x,[float(r[f'Tw_{key}_K']) for r in nodes],color=COLORS['tw'],ls=':',lw=1.7,label='Wall temperature')
        shell_lines=[]
        for j in range(n_shells):
            line,=taxis.plot(x,[float(r[f'Ts{j}_{key}_K']) for r in nodes],color=SOIL_COLORS[j%len(SOIL_COLORS)],lw=1.25,alpha=.9,label=f'Soil shell {j+1}')
            shell_lines.append(line)
        if legend_handles is None: legend_handles=[p,m,tg,tw]+shell_lines
        ax.set_title(title,loc='left',fontsize=10.5,fontweight='semibold')
        ax.set_ylabel('Pressure (MPa)',color=COLORS['p']); maxis.set_ylabel('Mass flow (kg/s)',color=COLORS['m']); taxis.set_ylabel('Temperature (K)',color=COLORS['tg'])
        ax.tick_params(axis='y',colors=COLORS['p']); maxis.tick_params(axis='y',colors=COLORS['m']); taxis.tick_params(axis='y',colors=COLORS['tg'])
        maxis.set_ylim(center-pad,center+pad); maxis.yaxis.set_major_locator(MaxNLocator(nbins=5)); maxis.yaxis.set_major_formatter(FormatStrFormatter('%.5f'))
        taxis.yaxis.set_major_locator(MaxNLocator(nbins=5)); taxis.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax.grid(True,color='#d9dee5',linewidth=.7,alpha=.85)
    z=[float(r['elevation_m']) for r in nodes]
    axes[4].fill_between(x,z,color=COLORS['z'],alpha=.2); axes[4].plot(x,z,color=COLORS['z'],lw=1.8)
    axes[4].set_title('Pipe elevation',loc='left',fontsize=10.5,fontweight='semibold'); axes[4].set_ylabel('Elevation (m)',color=COLORS['z']); axes[4].tick_params(axis='y',colors=COLORS['z']); axes[4].grid(True,color='#d9dee5',linewidth=.7,alpha=.85)
    axes[4].set_xlabel('Distance along pipe (km)'); axes[4].set_xlim(x[0],x[-1])
    fig.suptitle('Fully thermal pipe with radial soil shells: steady initializers',fontsize=14,fontweight='semibold',y=.997)
    fig.legend(legend_handles,[h.get_label() for h in legend_handles],loc='upper center',bbox_to_anchor=(.5,.979),ncol=min(6,len(legend_handles)),frameon=False,fontsize=9)
    fig.tight_layout(rect=(.04,.02,.75,.945),h_pad=1.0)
    fig.savefig(output,dpi=185,facecolor='white',bbox_inches='tight'); plt.close(fig)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('csv',type=Path,help='CSV from soil-shell profile exporter')
    parser.add_argument('-o','--output',type=Path,default=Path('soil_shell_initializer_profiles.png'))
    args=parser.parse_args(); plot_profiles(read_profiles(args.csv),args.output); print(f'Wrote {args.output}')

if __name__=='__main__': main()

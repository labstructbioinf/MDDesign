#! /usr/bin/env python
import argparse
import matplotlib
matplotlib.use('Agg')
import numpy as np
import seaborn as sns
import pandas as pd
import os
import sys
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import ListedColormap
from utils import get_fn, prHeader, prError, prOK, prWarning
parser = argparse.ArgumentParser(description='MDPlot')
parser.add_argument('-rmsd', nargs='+')
parser.add_argument('-secstruct', action='store_true')
parser.add_argument('-list')
parser.add_argument('-skip_first', type=int, default=0)
parser.add_argument('-skip_last', type=int, default=0)
parser.add_argument('-out', default='OUT.pdf')
parser.add_argument('-dt', type=float, default=2)
parser.add_argument('-rmsd_min', type=float)
parser.add_argument('-rmsd_max', type=float)
parser.add_argument('-percentile', type=float, default=100)
parser.add_argument('-format', default='png')
parser.add_argument('-dpi', default=300, type=float)
parser.add_argument('-linewidth', default=0.25, type=float)
parser.add_argument('-keep_images', action='store_true')
args = parser.parse_args()
header = ("""
       __   __        __  ___
 |\/| |  \ |__) |    /  \  |
 |  | |__/ |    |___ \__/  |
""")
prHeader(header)
sns.set(style='white')
sns.set_context(rc = {'patch.linewidth': 0.0})
rmsd_all = np.array
hist_all = []
if not any(os.access(os.path.join(path, 'convert'), os.X_OK) for path in os.environ["PATH"].split(os.pathsep)):
    prWarning("convert executable not found! Output will be separate images instead of merged PDF.")
    args.keep_images = True
if not os.path.isfile(args.list):
    prError("Model list does not exist.")
    sys.exit(1)
### Check if all files exist
missing_files = 0
with open(args.list) as f:
    for line in f:
        line = line.rstrip()
        s = line.split(';')
        ident = s[0]
        natoms = s[1]
        replicas = int(s[2])
        j = 1
        if args.rmsd:
            for suffix in args.rmsd:
                for i in range(1, replicas + 1):
                    if not os.path.isfile('./%s/%s_%s/%s_%s_%s.out' % (ident, ident, i, ident, i, suffix)):
                        prError('RMSD file for %s_%s (suffix %s) does not exist.' % (ident, i, suffix))
                        missing_files += 1
                    else:
                        z = np.genfromtxt('./%s/%s_%s/%s_%s_%s.out' % (ident, ident, i, ident, i, suffix))
                        if args.skip_first >= z.shape[0] or args.skip_last >= z.shape[0]:
                            prError(
                                "%s_%s (%s): Number of skipped frames (%s) is bigger than number of frames in the RMSD out file (%s)." % (
                                ident, i, suffix, (args.skip_last + args.skip_first), z.shape[0]))
                            missing_files += 1
        if args.secstruct:
             for i in range(1, replicas + 1):
                 if not os.path.isfile('./%s/%s_%s/%s_%s_secstruct.out' % (ident, ident, i, ident, i)):
                     prError('Secondary structure file for %s_%s does not exist.' % (ident, i))
                     missing_files += 1
                 else:
                     ss_df = pd.read_fwf('./%s/%s_%s/%s_%s_secstruct.out' % (ident, ident, i, ident, i))
                     if args.skip_first >= ss_df.shape[0] or args.skip_last >= ss_df.shape[0]:
                         prError(
                             "%s_%s (%s): Number of skipped frames (%s) is bigger than number of frames in the secondary structure out file (%s)." % (
                                 ident, i, suffix, (args.skip_last + args.skip_first), z.shape[0]))
                         missing_files += 1

if missing_files > 0:
    sys.exit(1)
else:
    prOK('All file checks passed.')
with open(args.list) as f:
    for line in f:
        line = line.rstrip()
        s = line.split(';')
        ident = s[0]
        natoms = s[1]
        replicas = int(s[2])
        j = 1
        if args.rmsd:
            for suffix in args.rmsd:
                k = np.array
                for i in range(1, replicas+1):
                    try:
                        z = np.genfromtxt('./%s/%s_%s/%s_%s_%s.out' % (ident, ident, i, ident, i, suffix),
                                          skip_header=args.skip_first + 1, skip_footer=args.skip_last + 1, usecols=(1,))

                    except IOError:

                        sys.exit(1)
                    if i == 1:
                        k = z
                    else:
                        k = np.concatenate((k, z))

                    #np.concatenate((hist_all, np.float(np.histogram(k, bins=100, density=True)[0].max())))
                   # print max
                    if j == 1:
                        rmsd_all = z
                    else:
                        rmsd_all = np.concatenate((rmsd_all, z))
                    j += 1
                axlim = sns.distplot(k, label="%s" % suffix)
                hist_all.append(float(axlim.get_ylim()[1]))
                plt.clf()

if args.rmsd_max:
    plot_ylim = args.rmsd_max
else:
    plot_ylim = np.percentile(rmsd_all, args.percentile) + 1
hist_ylim = max(hist_all)
with open(args.list) as f:
    for line in f:
        line = line.rstrip()
        s = line.split(';')
        ident = s[0]
        natoms = s[1]
        replicas = int(s[2])
        merged = np.array
        fig, ax = plt.subplots()
        G = GridSpec(2, 2)
        axes_hist = plt.subplot(G[0, 0])
        plt.suptitle('%s' % ident, fontsize=14, fontweight='bold')
        if args.rmsd:
            for suffix in args.rmsd:
                for i in range(1, replicas+1):
                    try:
                        z = np.genfromtxt('./%s/%s_%s/%s_%s_%s.out' % (ident, ident, i, ident, i, suffix),
                                          skip_header=args.skip_first + 1, skip_footer=args.skip_last, usecols=(1,))
                    except IOError:
                        pass
                    if i == 1:
                        merged = z
                    else:
                        merged = np.concatenate((merged, z))
                sns.distplot(merged, ax=axes_hist, label="%s" % suffix)
                axes_hist.set_xlim(0, plot_ylim)
                axes_hist.set_ylim(0, hist_ylim)
                axes_hist.set_xlabel('RMSD (A)', fontsize=8)
                leg = axes_hist.legend(fontsize=8)
            axes_hist.tick_params(labelsize=8, labelcolor="black")
            for legobj in leg.legendHandles:
                legobj.set_linewidth(2.0)
            axes_joint = plt.subplot(G[0, 1])
            for suffix in args.rmsd:
                for i in range(1, replicas+1):
                    z = np.genfromtxt('./%s/%s_%s/%s_%s_%s.out' % (ident, ident, i, ident, i, suffix),
                                      skip_header=args.skip_first + 1, skip_footer=args.skip_last, usecols=(1,))
                    if i == 1:
                        merged = z
                    else:
                        merged = np.concatenate((merged, z))
                axes_joint.plot(np.linspace(0, merged.size * (args.dt/float(1000)), num = merged.size), merged, label="%s" % suffix,
                                linewidth=0.25)
                axes_joint.set_xlim(0, merged.size * (args.dt/float(1000)))
                axes_joint.set_ylim(0, plot_ylim)
                axes_joint.set_ylabel('RMSD (A)', fontsize=8)
                axes_joint.set_xlabel('Time (ns)', fontsize=8)
                axes_joint.tick_params(labelsize=8, labelcolor="black")
                leg = axes_joint.legend(fontsize=8)
            for legobj in leg.legendHandles:
                legobj.set_linewidth(2.0)
        if args.secstruct:
            axes_ss = plt.subplot(G[1, 0])
            ss_df = pd.read_fwf('./%s/%s_%s/%s_%s_secstruct.out' % (ident, ident, i, ident, i))[args.skip_first:]
            for i in range(2, replicas+1):
                ss_df = pd.concat([ss_df, pd.read_fwf('./%s/%s_%s/%s_%s_secstruct.out' % (ident, ident, i, ident, i))[
                                          args.skip_first:]])
            data = np.rot90(np.array(ss_df)[:, 1:])
            palette = ["#FFFFFF", "#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#FFA500", "#2ecc71", "#34495e"]
            ax = sns.heatmap(data, ax=axes_ss, cbar_kws={"boundaries": np.linspace(-0.5, 7.5, 9)},
                             cmap=ListedColormap(palette), xticklabels=int(data.shape[1] / 5),
                             yticklabels=int(data.shape[0] / 7))
            colorbar = ax.collections[0].colorbar
            for axis in ['top', 'bottom', 'left', 'right']:
                axes_ss.spines[axis].set_visible(True)
                axes_ss.spines[axis].set_color('black')
            colorbar.set_ticks([0, 1, 2, 3, 4, 5, 6, 7])
            colorbar.set_ticklabels(['None', 'Para', 'Anti', '3-10', 'Alpha', 'Pi', 'Turn', 'Bend'])
            axes_ss.set_ylabel('Residue', fontsize=8)
            axes_ss.set_xlabel('Frames', fontsize=8)
            _, labels = plt.yticks()
            axes_ss.tick_params(labelsize=8, labelcolor="black")
        plt.tight_layout()
        plt.subplots_adjust(top=0.90)
        plt.savefig('%s.%s' % (ident, args.format), dpi=args.dpi)
        plt.close(fig)
        plt.clf()
convert_string = ''
with open(args.list) as f:
    for line in f:
        line = line.rstrip()
        s = line.split(';')
        ident = s[0]
        natoms = s[1]
        replicas = int(s[2])
        convert_string += '%s.%s ' % (ident, args.format)
os.system('convert %s %s' % (convert_string, args.out))
if not args.keep_images:
    os.system('rm *.%s' % args.format)
prOK("Plotting done!")

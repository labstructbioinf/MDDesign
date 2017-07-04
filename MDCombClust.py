#! /usr/bin/env python
import argparse
import sys
import subprocess
import os
from utils import get_fn, prHeader
script_path = os.path.dirname(__file__)
from utils import get_fn, prError, prOK, prWarning, prHeader
header = ("""
       __   __   __         __   __             __  ___
 |\/| |  \ /  ` /  \  |\/| |__) /  ` |    |  | /__`  |
 |  | |__/ \__, \__/  |  | |__) \__, |___ \__/ .__/  |
""")
prHeader(header)
parser = argparse.ArgumentParser(description='Analyzer')
parser.add_argument('-list')
parser.add_argument('-start_frame', type=int)
parser.add_argument('-last_frame', type=int)
parser.add_argument('-stride', type=int)
parser.add_argument('-step', type=int)
parser.add_argument('-nstruct', type=int)
parser.add_argument('-out_dir', default='COMBCLUST')
parser.add_argument('-ai', action='store_true')
args = parser.parse_args()
if not os.path.isfile(args.list):
    prError("Model list does not exist.")
    sys.exit(1)
def execute(command):
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    # Poll process for new output until finished
    while True:
        nextline = process.stdout.readline()
        if nextline == '' and process.poll() is not None:
            break
        if 'initial clusters' in nextline or '+' in nextline or 'Pair-wise matrix set up' in nextline or 'Memory us' in nextline:
            sys.stdout.write(nextline)
            sys.stdout.flush()

    output = process.communicate()[0]
    exitCode = process.returncode

    if (exitCode == 0):
        return output
errors = 0
if not any(os.access(os.path.join(path, 'cpptraj'), os.X_OK) for path in os.environ["PATH"].split(os.pathsep)):
    prError("cpptraj not found! Check amber installation!")
    errors += 1
if errors != 0:
    sys.exit(1)
try:
    os.mkdir(args.out_dir)
except OSError:
    if os.path.isdir(args.out_dir):
        if os.listdir(args.out_dir) == []:
            prWarning('Warning! Output directory exists and is not empty!')
print("Starting clustering with CPPTRAJ (this can take few hours)...")
f2 = open('combclust.in', 'w')
with open(args.list) as f:
    j = 0
    missing_files = 0
    mismatch_files = 0
    for line in f:
        file_sizes = set()
        j += 1
        line = line.rstrip()
        s = line.split(';')
        ident = s[0]
        natoms = s[1]
        replicas = int(s[2])
        f2.write('parm %s [%s]\n' % (get_fn(ident, 1, type='top'), ident,))
        for i in range(1, replicas + 1):
            if args.ai:
                if not os.path.isfile('%s/%s_%s/_%s/%s_%s_%s_ai.nc' % (ident, ident, i, args.step, ident, i, args.step)):
                    prError('Trajectory for model %s_%s does not exist.' % (ident, i))
                    missing_files += 1
                else:
                    file_sizes.add(
                        os.path.getsize('%s/%s_%s/_%s/%s_%s_%s_ai.nc' % (ident, ident, i, args.step, ident, i, args.step)))
            else:
                if not os.path.isfile('%s/%s_%s/_%s/%s_%s_%s.nc' % (ident, ident, i, args.step, ident, i, args.step)):
                    prError('Trajectory for model %s_%s does not exist.' % (ident, i))
                    missing_files += 1
                else:
                    file_sizes.add(
                        os.path.getsize('%s/%s_%s/_%s/%s_%s_%s.nc' % (ident, ident, i, args.step, ident, i, args.step)))
            if args.ai:
                f2.write('loadcrd %s %s %s %s parm [%s] name %s_trj\n' % (
                    get_fn(ident, i, step=args.step, type='traj_ai'), args.start_frame, args.last_frame, args.stride, ident, ident))
            else:
                f2.write('loadcrd %s %s %s %s parm [%s] name %s_trj\n' % (
                    get_fn(ident, i, step=args.step, type='traj'), args.start_frame, args.last_frame, args.stride,
                    ident, ident))
        f2.write('crdaction %s_trj strip !(@C,CA,O,N)\n' % (ident))
        f2.write('crdout %s_trj %s_trj_bb.nc\n' % (ident, ident))
        if j == 1:
            f2.write('parmwrite out bb.parm7 crdset %s_trj\n' % (ident))
        if len(file_sizes) > 1:
            prWarning('Mismatch in filesizes of trajectories for model %s' % ident)
            mismatch_files += 1
    f2.write('clear all\n')
if missing_files != 0:
    prError('Trajectories are missing. Rerun the Amber jobs for this step or exclude models above from the model list')
    sys.exit(1)
with open(args.list) as f:
    j = 0
    f2.write('parm bb.parm7\n')
    for line in f:
        line = line.rstrip()
        s = line.split(';')
        ident = s[0]
        natoms = s[1]
        replicas = int(s[2])
        f2.write('trajin %s_trj_bb.nc\n' % (ident))
    f2.write('rms first @C,CA,O,N\n')
    f2.write(
        'cluster out ./%s/cluster_out.out summary ./%s/cluster_summary.out info ./%s/cluster_info.out repout ./%s/rep_bb repfmt pdb averagelinkage clusters %s rms nofit @C,CA,O,N\n' % (
            args.out_dir, args.out_dir, args.out_dir, args.out_dir, args.nstruct))
    f2.close()
execute('cpptraj -i combclust.in')
os.system('rm combclust.in')
# Create a mapping from the frame number to the topology
with open('./%s/cluster_info.out' % args.out_dir) as f:
    for line in f:
        line = line.rstrip()
        if '#Representative frames:' in line:
            line = line.replace('#Representative frames: ', '')
            frames = line.split(' ')
mapping = {}
with open(args.list) as f:
    cnt = 1
    for line in f:
        line = line.rstrip()
        s = line.split(';')
        ident = s[0]
        natoms = s[1]
        replicas = int(s[2])
        top = '[%s]' % ident
        for i in range(1, replicas + 1):
            count = args.start_frame
            while count <= args.last_frame:
                count += args.stride
                mapping[cnt] = top
                cnt += 1
f2 = open('getfull.in', 'w')
with open(args.list) as f:
    for line in f:
        line = line.rstrip()
        s = line.split(';')
        ident = s[0]
        natoms = s[1]
        replicas = int(s[2])
        f2.write('parm ./%s/%s_1/%s_1.parm7 [%s]\n' % (ident, ident, ident, ident))
        for i in range(1, replicas + 1):
            if args.ai:
                f2.write('trajin %s %s %s %s parm [%s]\n' % (
                    get_fn(ident, i, step=args.step, type='traj_ai'), args.start_frame, args.last_frame, args.stride,
                    ident))
            else:
                f2.write('trajin %s %s %s %s parm [%s]\n' % (
                    get_fn(ident, i, step=args.step, type='traj'), args.start_frame, args.last_frame, args.stride,
                    ident))
    cnt = 1
    for frame in frames:
        f2.write('trajout ./%s/rep_allatom.c%s.pdb onlyframes %s parm %s\n' % (args.out_dir, cnt-1, frame, mapping[int(frame)]))
        cnt += 1
f2.close()
print("Getting full atom structures...")
execute('cpptraj -i getfull.in')
os.system('rm getfull.in')
os.system('rm bb.parm7 *_trj_bb.nc')
print("Done!")
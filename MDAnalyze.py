#! /usr/bin/env python
import argparse
import subprocess
import multiprocessing
import sys
import os
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
import shutil
parser = argparse.ArgumentParser(description='MDAnalyzer')
parser.add_argument('-job', help='Job type. Available are: autoimage, refstruct rmsd_bb, rmsd_sse, rmsd_custom, secstruct.', required = True)
parser.add_argument('-step', nargs='+', help='Simulation step to perform calculations on.', required=True)
parser.add_argument('-list', required=True, help='Model list')
parser.add_argument('-suffix', help='Suffix to identify rmsd calculation in plotting')
parser.add_argument('-mask', help='Amber mask for rmsd_custom calculation')
parser.add_argument('-nofit', action='store_true', help='Used with rmsd_custom. First backbone RMS is performed then for the range defined in -mask rmsd is performed without fitting.')
parser.add_argument('-ai', action='store_true')
parser.add_argument('-start_frame', default=1)
parser.add_argument('-end_frame', default=9999999)

from utils import get_fn, prError, prOK, prWarning, prHeader
args = parser.parse_args()


### Adapted from https://gist.github.com/aubricus/f91fb55dc6ba5557fbab06119420dd6a
def printProgressBar (iteration, total, prefix = '', suffix = ''):
    percent = ("{0:.0f}").format(100 * (iteration / float(total)))
    sys.stdout.write('\r%s %s%% %s' % (prefix, percent, suffix))
    sys.stdout.flush()
    # Print New Line on Complete
    if iteration == total:
        prOK("OK")
def os_cmd(cmd):
    subprocess.call([cmd], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #os.system(cmd)

def os_multiprocess(tasks):
    workers = multiprocessing.cpu_count()
    printProgressBar(0, len(tasks))
    stdout_queue = multiprocessing.Queue()
    pool = multiprocessing.Pool(
        processes=workers, initargs=[stdout_queue])
    for i, _ in enumerate(pool.imap(os_cmd, tasks), 1):
        printProgressBar(i, len(tasks))
    pool.close()
    pool.join()

header = ("""
       __                          __  ___
 |\/| |  \  /\  |\ |  /\  |    \ /  / |__
 |  | |__/ /~~\ | \| /~~\ |___  |  /_ |___
""")
prHeader(header)

errors = 0
if not any(os.access(os.path.join(path, 'cpptraj'), os.X_OK) for path in os.environ["PATH"].split(os.pathsep)):
    prError("cpptraj not found! Check amber installation!")
    errors += 1
if not any(os.access(os.path.join(path, 'tleap'), os.X_OK) for path in os.environ["PATH"].split(os.pathsep)):
    prError("tleap not found! Check amber installation!")
if not any(os.access(os.path.join(path, 'tleap'), os.X_OK) for path in os.environ["PATH"].split(os.pathsep)):
    prError("ante-MMPBSA.py not found! Check amber installation!")
    errors += 1
if errors != 0:
    sys.exit(1)
if not os.path.isfile(args.list):
    prError("Model list does not exist.")
    sys.exit(1)
### AUTOIMAGE ###

if args.job == 'autoimage':
    jobs_autoimage = []
    missing_files = 0
    mismatch_files = 0
    for step in args.step:
        # print step
        print ("Autoimage for step %s..." % (step))
        with open(args.list) as f:
            for line in f:
                file_sizes = set()

                line = line.rstrip()
                s = line.split(';')
                ident = s[0]
                replicas = int(s[2])

                for i in range(1, replicas + 1):
                    if not os.path.isfile('%s/%s_%s/_%s/%s_%s_%s.nc' % (ident, ident, i, step, ident, i, step)):
                        prError('Trajectory for model %s_%s does not exist.' % (ident, i))
                        missing_files += 1


                    else:
                        file_sizes.add(
                            os.path.getsize('%s/%s_%s/_%s/%s_%s_%s.nc' % (ident, ident, i, step, ident, i, step)))
                    f_input = open('./%s/%s_%s/%s_%s_ai_%s.in' % (ident, ident, i, ident, i, step), 'w')
                    inp = (
                        """
parm %s/%s_%s/%s_%s.parm7
trajin %s/%s_%s/_%s/%s_%s_%s.nc
autoimage
trajout %s/%s_%s/_%s/%s_%s_%s_ai.nc
""" % (ident, ident, i, ident, i, ident, ident, i, step, ident, i, step, ident, ident, i, step, ident, i, step)
                    )
                    f_input.write(inp)
                    f_input.close()
                    jobs_autoimage.append('cpptraj -i ./%s/%s_%s/%s_%s_ai_%s.in' % (ident, ident, i, ident, i, step))
                if len(file_sizes) > 1:
                    prWarning('Mismatch in filesizes of trajectories for model %s' % ident)
                    mismatch_files += 1

    if missing_files != 0:
        prError('Trajectories are missing. Rerun the Amber jobs for this step or exclude models above from the model list')
        sys.exit(1)
    if mismatch_files != 0:
        prWarning('Some files may be corrupted.')
    os_multiprocess(jobs_autoimage)

### REFSTRUCT ###

elif args.job == 'refstruct':
    jobs_refstruct = []
    missing_files = 0
    print ("Preparing reference structure (restart from step %s)..." % (args.step[0]))
    with open(args.list) as f:
        for line in f:
            line = line.rstrip()
            s = line.split(';')
            ident = s[0]
            natoms = s[1]
            replicas = int(s[2])
            for i in range(1, replicas + 1):
                subprocess.call(['rm ./%s/%s_%s/%s_%s_ref.rst7' % (ident, ident, i, ident, i)], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                f_input = open('./%s/%s_%s/%s_%s_ref.in' % (ident, ident, i, ident, args.step[0]), 'w')
                inp = (
                          """
parm %s/%s_%s/%s_%s_i_s.parm7
trajin %s/%s_%s/_%s/%s_%s_%s.rst7
strip :WAT,Na+,Cl-
trajout ./%s/%s_%s/%s_%s_ref.rst7
""") % (ident, ident, i, ident, i, ident, ident, i, args.step[0], ident, i, args.step[0], ident, ident, i, ident, i)
                f_input.write(inp)
                f_input.close()
                jobs_refstruct.append("cpptraj -i ./%s/%s_%s/%s_%s_ref.in" % (ident, ident, i, ident, args.step[0]))
    os_multiprocess(jobs_refstruct)
    with open(args.list) as f:
        for line in f:
            line = line.rstrip()
            s = line.split(';')
            ident = s[0]
            natoms = s[1]
            replicas = int(s[2])
            for i in range(1, replicas + 1):
                if not os.path.isfile('./%s/%s_%s/%s_%s_ref.rst7' % (ident, ident, i, ident, i)):
                    prError('Error: reference structure for model %s_%s was not generated' % (ident, i))
                    missing_files += 1
    if missing_files != 0:
        prError('Errors in reference structure generation occured.')
        prError('Did you propagate files for first step? e.g. MDAnalyze -list MODEL_LIST -copy -step 1')
        sys.exit(1)
    else:
        prOK("All files generated successfully!")

### RMSD ###

elif args.job == 'rmsd_bb':
    missing_files = 0
    file_sizes = set()
    jobs_rmsd = []
    step_string = ""
    for step in args.step:
        step_string += "%s " % step
    print ("Backbone RMSD analysis for steps: %s" % (step_string))
    with open(args.list) as f:
        for line in f:
            line = line.rstrip()
            s = line.split(';')
            ident = s[0]
            natoms = s[1]
            replicas = int(s[2])
            for i in range(1, replicas + 1):
                subprocess.call(['rm %s/%s_%s/%s_%s_%s.out' % ( ident, ident, i, ident, i, args.suffix)], shell=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                f_input = open('./%s/%s_%s/%s_%s_%s.in' % (ident, ident, i, ident, i, args.suffix), 'w')
                input = "parm %s/%s_%s/%s_%s.parm7\n" % (ident, ident, i, ident, i)
                input += "reference %s/%s_%s/%s_%s_ref.rst7\n" % (ident, ident, i, ident, i)
                k = 0
                for step in args.step:
                    if args.ai:
                        input += "trajin %s/%s_%s/_%s/%s_%s_%s_ai.nc\n" % (ident, ident, i, step, ident, i, step)
                    else:
                        input += "trajin %s/%s_%s/_%s/%s_%s_%s.nc\n" % (ident, ident, i, step, ident, i, step)
                input += "rms reference @C,CA,O,N out %s/%s_%s/%s_%s_%s.out \n" % (
                    ident, ident, i, ident, i, args.suffix)
                f_input.write(input)
                f_input.close()
                jobs_rmsd.append("cpptraj -i ./%s/%s_%s/%s_%s_%s.in" % (ident, ident, i, ident, i, args.suffix))
    os_multiprocess(jobs_rmsd)
    with open(args.list) as f:
        for line in f:
            line = line.rstrip()
            s = line.split(';')
            ident = s[0]
            natoms = s[1]
            replicas = int(s[2])
            for i in range(1, replicas + 1):
                if not os.path.isfile('%s/%s_%s/%s_%s_%s.out' % (
                    ident, ident, i, ident, i, args.suffix)):
                    prError('Error: RMSD output for model %s_%s was not generated' % (ident, i))
                    missing_files += 1
                else:
                    file_sizes.add(
                        os.path.getsize('%s/%s_%s/%s_%s_%s.out' % (
                    ident, ident, i, ident, i, args.suffix)))
    if missing_files != 0:
        prError('Errors in RMSD calculation occured.')
        prError('Please autoimage trajectory, calculate reference structure or rerun this step in Amber.')
        sys.exit(1)
    elif len(file_sizes) != 1:
        prWarning('Mismatch in generated file sizes. Files may be corrupt.')
    else:
        prOK("All files generated successfully!")

elif args.job == 'rmsd_heavy':
    missing_files = 0
    file_sizes = set()
    jobs_rmsd = []
    step_string = ""
    for step in args.step:
        step_string += "%s " % step
    print ("Heavy atoms RMSD analysis for steps: %s" % (step_string))
    with open(args.list) as f:
        for line in f:
            line = line.rstrip()
            s = line.split(';')
            ident = s[0]
            natoms = s[1]
            replicas = int(s[2])
            for i in range(1, replicas + 1):
                subprocess.call(['rm %s/%s_%s/%s_%s_%s.out' % (ident, ident, i, ident, i, args.suffix)], shell=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                f_input = open('./%s/%s_%s/%s_%s_%s.in' % (ident, ident, i, ident, i, args.suffix), 'w')
                input = "parm %s/%s_%s/%s_%s.parm7\n" % (ident, ident, i, ident, i)
                input += "reference %s/%s_%s/%s_%s_ref.rst7\n" % (ident, ident, i, ident, i)
                k = 0
                for step in args.step:
                    if args.ai:
                        input += "trajin %s/%s_%s/_%s/%s_%s_%s_ai.nc\n" % (ident, ident, i, step, ident, i, step)
                    else:
                        input += "trajin %s/%s_%s/_%s/%s_%s_%s.nc\n" % (ident, ident, i, step, ident, i, step)
                input += "rms reference !@H= out %s/%s_%s/%s_%s_%s.out \n" % (ident, ident, i, ident, i, args.suffix)
                f_input.write(input)
                f_input.close()
                jobs_rmsd.append("cpptraj -i ./%s/%s_%s/%s_%s_%s.in" % (ident, ident, i, ident, i, args.suffix))
    os_multiprocess(jobs_rmsd)
    with open(args.list) as f:
        for line in f:
            line = line.rstrip()
            s = line.split(';')
            ident = s[0]
            natoms = s[1]
            replicas = int(s[2])
            for i in range(1, replicas + 1):
                if not os.path.isfile('%s/%s_%s/%s_%s_%s.out' % (
                    ident, ident, i, ident, i, args.suffix)):
                    prError('Error: RMSD output for model %s_%s was not generated' % (ident, i))
                    missing_files += 1
                else:
                    file_sizes.add(
                        os.path.getsize('%s/%s_%s/%s_%s_%s.out' % (
                    ident, ident, i, ident, i, args.suffix)))
    if missing_files != 0:
        prError('Errors in RMSD calculation occured.')
        prError('Please autoimage trajectory, calculate reference structure or rerun this step in Amber.')
        sys.exit(1)
    elif len(file_sizes) != 1:
        prWarning('Mismatch in generated file sizes. Files may be corrupt.')
    else:
        prOK("All files generated successfully!")

elif args.job == 'rmsd_sse':
    jobs_rmsd = []
    missing_files = 0
    file_sizes = set()
    step_string = ""
    for step in args.step:
        step_string += "%s " % step
    print ("RMSD over the secondary structure elements analysis for steps: %s" % (step_string))
    with open(args.list) as f:
        for line in f:
            line = line.rstrip()
            s = line.split(';')
            ident = s[0]
            replicas = int(s[2])
            for i in range(1, replicas + 1):
                subprocess.call(['rm %s/%s_%s/%s_%s_%s.out' % (ident, ident, i, ident, i, args.suffix)], shell=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                f_input = open('./%s/%s_%s/%s_%s_%s.in' % (ident, ident, i, ident, i, args.suffix), 'w')
                input = "parm %s/%s_%s/%s_%s.parm7\n" % (ident, ident, i, ident, i)
                input += "reference %s/%s_%s/%s_%s_ref.rst7\n" % (ident, ident, i, ident, i)
                k = 0
                for step in args.step:
                    if args.ai:
                        input += "trajin %s/%s_%s/_%s/%s_%s_%s_ai.nc\n" % (ident, ident, i, step, ident, i, step)
                    else:
                        input += "trajin %s/%s_%s/_%s/%s_%s_%s_ai.nc\n" % (ident, ident, i, step, ident, i, step)
                cmd = "ambpdb -bres -p %s/%s_%s/%s_%s.parm7 -c %s/%s_%s/%s_%s_ref.rst7  > temp.pdb" % (
                    ident, ident, i, ident, i, ident, ident, i, ident, i)
                subprocess.call([cmd], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                dssp_dict, keys = dssp_dict_from_pdb_file("temp.pdb")
                acceptable_ss = ['E', 'H']
                sse_residues = ""
                for key in keys:
                    residue = key[1][1]
                    ss = dssp_dict[key][1]
                    if ss in acceptable_ss:
                        sse_residues += ":"
                        sse_residues += str(residue)
                        sse_residues += ","

                sse_residues = sse_residues[:-1]
                input += "rms reference %s&!@H= out %s/%s_%s/%s_%s_%s.out \n" % (
                    sse_residues, ident, ident, i, ident, i, args.suffix)
                f_input.write(input)
                f_input.close()
                jobs_rmsd.append("cpptraj -i ./%s/%s_%s/%s_%s_%s.in" % (ident, ident, i, ident, i, args.suffix))
                cmd = "rm temp.pdb"
                subprocess.call([cmd], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    os_multiprocess(jobs_rmsd)
    with open(args.list) as f:
        for line in f:
            line = line.rstrip()
            s = line.split(';')
            ident = s[0]
            natoms = s[1]
            replicas = int(s[2])
            for i in range(1, replicas + 1):
                if not os.path.isfile('%s/%s_%s/%s_%s_%s.out' % (
                    ident, ident, i, ident, i, args.suffix)):
                    prError('Error: RMSD output for model %s_%s was not generated' % (ident, i))
                    missing_files += 1
                else:
                    file_sizes.add(
                        os.path.getsize('%s/%s_%s/%s_%s_%s.out' % (
                    ident, ident, i, ident, i, args.suffix)))
    if missing_files != 0:
        prError('Errors in RMSD calculation occured.')
        prError('Please autoimage trajectory, calculate reference structure or rerun this step in Amber.')
        sys.exit(1)
    elif len(file_sizes) != 1:
        prWarning('Mismatch in generated file sizes. Files may be corrupt.')
    else:
        prOK("All files generated successfully!")
elif args.job == 'rmsd_custom':
    jobs_rmsd = []
    file_sizes = set()
    missing_files = 0
    step_string = ""
    for step in args.step:
        step_string += "%s " % step
    print ("Custom RMSD analysis for steps: %s" % (step_string))
    if args.nofit:
        with open(args.list) as f:
            for line in f:
                line = line.rstrip()
                s = line.split(';')
                ident = s[0]
                natoms = s[1]
                replicas = int(s[2])
                for i in range(1, replicas + 1):
                    subprocess.call(['rm %s/%s_%s/%s_%s_%s.out' % (ident, ident, i, ident, i, args.suffix)], shell=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    f_input = open('./%s/%s_%s/%s_%s_%s.in' % (ident, ident, i, ident, i, args.suffix), 'w')
                    input = "parm %s/%s_%s/%s_%s.parm7\n" % (ident, ident, i, ident, i)
                    input += "reference %s/%s_%s/%s_%s_ref.rst7\n" % (ident, ident, i, ident, i)
                    k = 0
                    for step in args.step:
                        if args.ai:
                            input += "trajin %s/%s_%s/_%s/%s_%s_%s_ai.nc\n" % (ident, ident, i, step, ident, i, step)
                        else:
                            input += "trajin %s/%s_%s/_%s/%s_%s_%s.nc\n" % (ident, ident, i, step, ident, i, step)
                    input += "rms reference @C,CA,O,N \n"
                    input += "rmsd %s out %s/%s_%s/%s_%s_%s.out nofit\n" % (
                        args.mask, ident, ident, i, ident, i, args.suffix)
                    f_input.write(input)
                    f_input.close()
                    jobs_rmsd.append("cpptraj -i ./%s/%s_%s/%s_%s_%s.in" % (ident, ident, i, ident, i, args.suffix))
    else:
        with open(args.list) as f:
            for line in f:
                line = line.rstrip()
                s = line.split(';')
                ident = s[0]
                natoms = s[1]
                replicas = int(s[2])
                for i in range(1, replicas + 1):
                    subprocess.call(['rm %s/%s_%s/%s_%s_%s.out' % (ident, ident, i, ident, i, args.suffix)], shell=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    f_input = open('./%s/%s_%s/%s_%s_%s.in' % (ident, ident, i, ident, i, args.suffix), 'w')
                    input = "parm %s/%s_%s/%s_%s.parm7\n" % (ident, ident, i, ident, i)
                    input += "reference %s/%s_%s/%s_%s_ref.rst7\n" % (ident, ident, i, ident, i)
                    k = 0
                    for step in args.step:
                        if args.ai:
                            input += "trajin %s/%s_%s/_%s/%s_%s_%s_ai.nc\n" % (ident, ident, i, step, ident, i, step)
                        else:
                            input += "trajin %s/%s_%s/_%s/%s_%s_%s.nc\n" % (ident, ident, i, step, ident, i, step)
                    input += "rms reference %s out %s/%s_%s/%s_%s_%s.out \n" % (
                        args.mask, ident, ident, i, ident, i, args.suffix)
                    f_input.write(input)
                    f_input.close()
                    jobs_rmsd.append("cpptraj -i ./%s/%s_%s/%s_%s_%s.in" % (ident, ident, i, ident, i, args.suffix))
    os_multiprocess(jobs_rmsd)
    with open(args.list) as f:
        for line in f:
            line = line.rstrip()
            s = line.split(';')
            ident = s[0]
            natoms = s[1]
            replicas = int(s[2])
            for i in range(1, replicas + 1):
                if not os.path.isfile('%s/%s_%s/%s_%s_%s.out' % (
                        ident, ident, i, ident, i, args.suffix)):
                    prError('Error: RMSD output for model %s_%s was not generated' % (ident, i))
                    missing_files += 1
                else:
                    file_sizes.add(
                        os.path.getsize('%s/%s_%s/%s_%s_%s.out' % (
                            ident, ident, i, ident, i, args.suffix)))
    if missing_files != 0:
        prError('Errors in RMSD calculation occured.')
        prError('Please autoimage trajectory, calculate reference structure or rerun this step in Amber.')
        sys.exit(1)
    elif len(file_sizes) != 1:
        prWarning('Mismatch in generated file sizes. Files may be corrupt.')
    else:
        prOK("All files generated successfully!")
### SECSTRUCT ###

elif args.job == 'secstruct':
    jobs_rmsd = []
    file_sizes = set()
    step_string = ""
    missing_files = 0
    for step in args.step:
        step_string += "%s " % step
    print ("Secondary structure analysis for steps: %s" % (step_string))
    with open(args.list) as f:
        for line in f:
            line = line.rstrip()
            s = line.split(';')
            ident = s[0]
            natoms = s[1]
            replicas = int(s[2])
            for i in range(1, replicas + 1):
                subprocess.call(['rm %s/%s_%s/%s_%s_secstruct.out' % (ident, ident, i, ident, i)], shell=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                f_input = open('./%s/%s_%s/%s_%s_secstruct.in' % (ident, ident, i, ident, i), 'w')
                input = "parm %s/%s_%s/%s_%s.parm7\n" % (ident, ident, i, ident, i)
                k = 0
                for step in args.step:
                    if args.ai:
                        input += "trajin %s/%s_%s/_%s/%s_%s_%s_ai.nc\n" % (ident, ident, i, step, ident, i, step)
                    else:
                        input += "trajin %s/%s_%s/_%s/%s_%s_%s.nc\n" % (ident, ident, i, step, ident, i, step)
                input += "secstruct out %s/%s_%s/%s_%s_secstruct.out \n" % (ident, ident, i, ident, i)
                f_input.write(input)
                f_input.close()
                jobs_rmsd.append("cpptraj -i ./%s/%s_%s/%s_%s_secstruct.in" % (ident, ident, i, ident, i))
    os_multiprocess(jobs_rmsd)
    with open(args.list) as f:
        for line in f:
            line = line.rstrip()
            s = line.split(';')
            ident = s[0]
            natoms = s[1]
            replicas = int(s[2])
            for i in range(1, replicas + 1):
                if not os.path.isfile('%s/%s_%s/%s_%s_secstruct.out' % (
                        ident, ident, i, ident, i)):
                    prError('Error: Secondary structure output for model %s_%s was not generated' % (ident, i))
                    missing_files += 1
                else:
                    file_sizes.add(
                        os.path.getsize('%s/%s_%s/%s_%s_secstruct.out' % (
                        ident, ident, i, ident, i)))
    if missing_files != 0:
        prError('Errors in secondary structure calculation occured.')
        prError('Please autoimage trajectory or rerun this step in Amber.')
        sys.exit(1)
    if len(file_sizes) != 1:
        prWarning('Mismatch in generated file sizes. Files may be corrupt.')
    else:
        prOK("All files generated successfully!")

else:
    prError("No or invalid job specified!")
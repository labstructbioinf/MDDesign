#! /usr/bin/env python
import os
import argparse
import subprocess
import multiprocessing
import re
import sys
from utils import get_fn, prError, prOK, prWarning, prHeader
script_path = os.path.dirname(__file__)


parser = argparse.ArgumentParser(description='MDPrep')
parser.add_argument('-rep', type=int, help='Number of individual replicas of simulations.')
parser.add_argument('-start', action='store_true', help='Used to initialize a directory tree a prepare first input.')
parser.add_argument('-i', default="", help="Amber input file for the current step.")
parser.add_argument('-list', default="", required=True, help="List with model information. If used with -start the list will be created, otherwise list will be parsed.")
parser.add_argument('-step', help="Current step in the simulation protocol. If used with -start '1' is used by default")
parser.add_argument('-one', action='store_true', help="Prepare inputs only for one replica. Useful for minimization.")
parser.add_argument('-copy', action='store_true', help="Propagate files after run with only one replica.")
parser.add_argument('-cmd', default='pmemd.cuda', help="Amber executable, default: 'pmemd.cuda'")
parser.add_argument('-nox', action='store_true', help="Do not output trajectory.")
parser.add_argument('-cluster_header', help="Header of the batch job used in cluster queue manager.")
parser.add_argument('-cluster_add_cmd', default='sbatch', help="Command adding batch job to the queue on the cluster. Default - 'sbatch'")
parser.add_argument('-ref', action='store_true', help="Use previous restart as reference structure. Useful when restraints are used.")
parser.add_argument('-sol_wall_dist', type=float, default=10.0, help="Solute-wall distance. Used in tleap for molecule solvation")
parser.add_argument('-use_top', action='store_true', help="Use provided (solvated!) topologies and coordinates instead of PDBS. In this case those should be provided as ID_i_s.parm7 and ID_i_s.rst7.")
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


def check_files(prefix, extensions):
    k = []
    ext = []
    c = 0
    for extension in extensions:
        k.append(os.path.isfile('%s%s' % (prefix, extension)))
        if os.path.isfile('%s%s' % (prefix, extension)) == True:
            ext.append([prefix, extension])
    if any(z == True for z in k):
        if len(ext) == 1:
            return ext
        else:
            return 2  #### Raise exception
    else:
        return 1


######################
#### Initialize ######

header = ("""
       __   __   __   ___  __
 |\/| |  \ |__) |__) |__  |__)
 |  | |__/ |    |  \ |___ |
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
if args.start:
    if not args.i:
        prError("Please provide Amber input. Usually first step is minimization.")
        sys.exit(1)
    if not os.path.isfile(args.i):
        prError("Amber input file does not exist.")
        sys.exit(1)
    print("Starting new directory tree.\n")
    # if not os.path.isfile(args.i):
    # print "Input file does not exist. Exiting..."
    # sys.exit(1)
    args_mdprep = []
    args_mdgen = []
    args_antemm = []
    args_rename = []
    args_move = []
    # Start from provided topologies
    if args.use_top:
        top_exts = (".parm", ".parm7", ".top", ".top7")
        crd_exts = (".crd", ".crd7", ".rst7", ".rst")
        tops = [f for f in os.listdir('./') if f.endswith(top_exts)]
        tops_proper = [os.path.splitext(os.path.basename(top))[0] for top in tops if
                       os.path.splitext(os.path.basename(top))[0].endswith("_i_s")]
        rsts = [f for f in os.listdir('./') if (f.endswith(crd_exts) and f.startswith(tuple(tops_proper)))]
        models = [os.path.splitext(os.path.basename(rst))[0].replace('_i_s', '') for rst in rsts if
                  os.path.splitext(os.path.basename(rst))[0].endswith("_i_s")]

        if len(tops) != 0:
            print("Found %s topologies and coordinates:" % len(tops))
            for model in models:
                # print model
                a = check_files("%s_i_s" % model, crd_exts)
                b = check_files("%s_i_s" % model, top_exts)
                if (a != 0 and a != 1) and (b != 0 and b != 1):
                    print(model)
                else:
                    models.remove(model)
            print(", ".join(tops))
            print ("")
        else:
            prError("Error: No topology and coordinate files found.")
            sys.exit(1)
        if not args.rep:
            prError("Error: Number of replicas not specified. Please use -rep N")
            sys.exit(1)
        try:
            os.mkdir('./_SH/')
        except OSError:
            pass
        for model in models:
            try:
                os.mkdir(model)
            except OSError:
                pass
            for i in range(1, args.rep + 1):

                try:
                    os.mkdir('./%s/%s_%s/' % (model, model, i))
                except OSError:
                    pass
                try:
                    os.mkdir('./%s/%s_%s/_1' % (model, model, i))
                except OSError:
                    pass
                try:
                    os.mkdir('./%s/%s_%s/_SH' % (model, model, i))
                except OSError:
                    pass
            args_antemm.append("ante-MMPBSA.py -p %s -s :WAT,:Na+,:Cl- -c %s" % (
                get_fn(model, 1, type='solv_top'), get_fn(model, 1, type='top')))
            args_mdgen.append(
                "python %s/gen.py -start -wd ./%s/%s_1/ -top %s_i_s%s -crd %s_i_s%s -id %s_1 -i %s -nox -exe '%s'" % (
                    script_path, model, model, model, check_files("%s_i_s" % model, top_exts)[0][1], model,
                    check_files("%s_i_s" % model, crd_exts)[0][1], model, args.i, args.cmd))
        print ("Preparing folder tree and minimization inputs...")
        os_multiprocess(args_mdgen)
        print ("Preparing stripped topologies for visualization...")
        os_multiprocess(args_antemm)
        pdbs = models

    # Start from PDBS
    else:
        pdbs = [f for f in os.listdir('./') if f.endswith('.pdb')]  ## Get ald PDB's in current directory
        if len(pdbs) != 0:
            print("Found %s PDB files:" % len(pdbs))
            print(", ".join(pdbs))
            print ("")
        else:
            prError("Error: No PDB files found in the working directory.")
            sys.exit(1)
        if not args.rep:
            prError("Error: Number of replicas not specified. Please use -rep N")
            sys.exit(1)
        try:
            os.mkdir('./_SH/')
        except OSError:
            pass
        for filename in pdbs:
            # Get model ID from the pdb filename (strip the .pdb ext)
            ident = os.path.splitext(os.path.basename(filename))[0]
            try:
                os.mkdir(ident)
            except OSError:
                pass
            args_mdprep.append(
                'python %s/prep.py -pdb %s -parm %s.parm7 -crd %s.rst7 -d %s -ions' % (
                    script_path, filename, ident, ident, args.sol_wall_dist))
            args_mdgen.append(
                "python %s/gen.py -start -wd ./%s/%s_1/ -top %s.parm7 -crd %s.rst7 -id %s_1 -i %s -nox -exe '%s'" % (
                    script_path, ident, ident, ident, ident, ident, args.i, args.cmd))
            args_antemm.append(
                "ante-MMPBSA.py -p %s -s :WAT,:Na+,:Cl- -c %s" % (
                    get_fn(ident, 1, type='solv_top'), get_fn(ident, 1, type='top')))
            for i in range(1, args.rep + 1):
                try:
                    os.mkdir('./%s/%s_%s/' % (ident, ident, i))
                except OSError:
                    pass
                try:
                    os.mkdir('./%s/%s_%s/_1' % (ident, ident, i))
                except OSError:
                    pass
                try:
                    os.mkdir('./%s/%s_%s/_SH' % (ident, ident, i))
                except OSError:
                    pass
        print("Preparing topologies and folder tree...")
        os_multiprocess(args_mdprep)
        print("Preparing folder tree and minimization inputs...")
        os_multiprocess(args_mdgen)
        ## Now check if all topologies were correctly generated
        error_count = 0
        for filename in pdbs:
            ident = os.path.splitext(os.path.basename(filename))[0]
            if os.path.isfile(get_fn(ident, 1, type='solv_top')) and os.path.isfile(get_fn(ident, 1, type='rst')):
                pass
            else:
                error_count += 1
                prWarning("Topology could not be created automatically for %s." % (ident))

        if error_count != 0:
            for filename in pdbs:
                ident = os.path.splitext(os.path.basename(filename))[0]
                subprocess.call(['rm -r %s' % (ident)], shell=True, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
            subprocess.call(['rm -r *.rst7 *.parm7 *leap.in *leap.log _SH'], shell=True, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
            prError(
                "\nTopology generation failed for at least one model.\nTry to remove any non-standard residues, waters etc. and rerun script or prepare topologies manually.\n")
            sys.exit(1)

        print("Preparing stripped topologies for visualization...")
        os_multiprocess(args_antemm)
    fml = open(args.list, 'w')
    fml.close()
    for filename in pdbs:
        ident = os.path.splitext(os.path.basename(filename))[0]
        subprocess.call(['rm *.rst7 *.parm7 *leap.in *leap.log'], shell=True, stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)
        for i in range(2, args.rep + 1):
            args_rename.append('cp %s %s' % (get_fn(ident, 1, type='solv_top'), get_fn(ident, i, type='solv_top')))
            args_rename.append('cp %s %s' % (get_fn(ident, 1, type='top'), get_fn(ident, i, type='top')))
            args_rename.append('cp %s %s' % (get_fn(ident, 1, type='rst'), get_fn(ident, i, type='rst')))
            args_rename.append(
                'cp %s %s' % (get_fn(ident, 1, step=1, type='inp'), get_fn(ident, i, step=1, type='inp')))
            args_rename.append(
                'cp %s %s' % (get_fn(ident, 1, step=1, type='sh'), get_fn(ident, i, step=1, type='sh')))
            args_rename.append(
                'cp %s %s' % (get_fn(ident, 1, type='mdg'), get_fn(ident, i, type='mdg')))
        with open("%s" % (get_fn(ident, 1, type='top')), 'r+') as f:
            lines = f.readlines()
            for i in range(0, len(lines)):
                if "FLAG POINTERS" in lines[i]:
                    s = re.split(r'\s+', lines[i + 2])
                    ntwprt = s[1]

        towrite = '%s;%s;%s\n' % (ident, ntwprt, args.rep)
        fml = open(args.list, 'a')
        fml.write(towrite)
        fml.close()
    if args.rep > 1:
        print ("Copying necessary files...")
        os_multiprocess(args_rename)
    with open(args.list) as f:
        f_batch = open('STEP_1.sh', 'w')
        for line in f:
            line = line.rstrip()
            s = line.split(';')
            ident = s[0]
            natoms = s[1]
            replicas = int(s[2])
            if args.one:
                replicas = 1
            for i in range(1, replicas + 1):
                sed_cmd = 'sed -i -- \'s/%s_1/%s_%s/g\' ./%s/%s_%s/%s_%s.mdg' % (
                    ident, ident, i, ident, ident, i, ident, i)
                subprocess.call([sed_cmd], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                sed_cmd = 'sed -i -- \'s/%s_1/%s_%s/g\' ./%s/%s_%s/_SH/1.sh' % (
                    ident,  ident, i, ident, ident, i,)
                subprocess.call([sed_cmd], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                fn = './_SH/%s_%s_1.sh' % (ident, i)
                jobid = '%s_%s' % (i, ident)
                inp = ''
                if args.cluster_header:
                    inp += open(args.cluster_header, 'r').read()
                inp += "cd %s/%s_%s\n" % (ident, ident, i)
                f2 = open('%s/%s_%s/_SH/1.sh' % (ident, ident, i), 'r')
                command = f2.readlines()[0].rstrip()
                f2.close()
                inp += "%s\n" % command
                f2 = open(fn, 'w')
                f2.write(inp)
                f2.close()
                if args.cluster_header:
                    f_batch.write(
                        '%s ./_SH/%s_%s_1.sh\nrm ./_SH/%s_%s_1.sh\n' % (args.cluster_add_cmd, ident, i, ident, i))
                else:
                    subprocess.call(['chmod +rx ./_SH/%s_%s_1.sh' % (ident, i)], shell=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
                    f_batch.write(
                        'echo MODEL %s_%s\n./_SH/%s_%s_1.sh\n' % (ident, i, ident, i))
        f_batch.close()
        subprocess.call(['chmod +rx STEP_1.sh'], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
elif args.copy:
    args_copy = []
    args_rename = []
    with open(args.list) as f:
        for line in f:
            line = line.rstrip()
            s = line.split(';')
            ident = s[0]
            natoms = s[1]
            replicas = int(s[2])
            for i in range(2, replicas + 1):
                try:
                    os.mkdir('%s/%s_%s/_%s' % (ident, ident, i, args.step))
                except OSError:
                    pass
                args_copy.append(
                    'cp %s %s' % (
                    get_fn(ident, 1, step=args.step, type='rst'), get_fn(ident, i, step=args.step, type='rst')))
                os.system(
                    'cp %s/%s_1/_SH/%s.sh %s/%s_%s/_SH/%s.sh' % (ident, ident, args.step, ident, ident, i, args.step))
                f = open('./%s/%s_%s/_SH/%s.sh' % (ident, ident, i, args.step), 'r')
                command = f.readlines()[0].rstrip()
                f.close()
                command = command.replace('_1.', '_%s.' % (args.step))
                command = command.replace('_1_', '_%s_' % (i))
                f = open('./%s/%s_%s/_SH/%s.sh' % (ident, ident, i, args.step), 'w')
                f.write(command)
                f.close()
    print ("Propagating files for step %s..." % (args.step))
    os_multiprocess(args_copy)
else:
    if not args.i:
        prError("Please provide Amber input (-i input.in).")
        sys.exit(1)
    if not os.path.isfile(args.i):
        prError("Amber input file does not exist.")
        sys.exit(1)
    if not os.path.isfile(args.list):
        prError("Model list does not exist.")
        sys.exit(1)
    if not args.step:
        prError("Please provide step number (-step N).")
        sys.exit(1)
    args_mdgen = []
    with open(args.list) as f:
        for line in f:
            line = line.rstrip()
            s = line.split(';')
            ident = s[0]
            natoms = s[1]
            replicas = int(s[2])
            if args.one:
                replicas = 1
            for i in range(1, replicas + 1):
                mdgen_cmd = "python %s/gen.py -wd %s -s %s -id %s_%s -i %s -exe '%s' -ntwprt %s " % (
                    script_path, get_fn(ident, i), args.step, ident, i, args.i, args.cmd, natoms)
                if args.nox:
                    mdgen_cmd += "-nox "
                if args.ref:
                    mdgen_cmd += "-ref "
                args_mdgen.append(mdgen_cmd)
        print ("Preparing inputs for step %s..." % (args.step))
    os_multiprocess(args_mdgen)
    with open(args.list) as f:
        f_batch = open('STEP_%s.sh' % (args.step), 'w')
        for line in f:
            line = line.rstrip()
            s = line.split(';')
            ident = s[0]
            natoms = s[1]
            replicas = int(s[2])
            if args.one:
                replicas = 1
            for i in range(1, replicas + 1):
                fn = './_SH/%s_%s_%s.sh' % (ident, i, args.step)
                jobid = '%s_%s' % (i, ident)
                inp = ''
                if args.cluster_header:
                    inp += open(args.cluster_header, 'r').read()
                inp += "cd %s/%s_%s\n" % (ident, ident, i)
                f2 = open('%s/%s_%s/_SH/%s.sh' % (ident, ident, i, args.step), 'r')
                command = f2.readlines()[0].rstrip()
                f2.close()
                inp += "%s\n" % (command)
                f2 = open(fn, 'w')
                f2.write(inp)
                f2.close()
                if args.cluster_header:
                    f_batch.write(
                        '%s ./_SH/%s_%s_%s.sh\nrm ./_SH/%s_%s_%s.sh\n' % (args.cluster_add_cmd, ident, i, args.step, ident, i, args.step))
                else:
                    subprocess.call(['chmod +rx ./_SH/%s_%s_%s.sh' % (ident, i, args.step)], shell=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
                    f_batch.write(
                        'echo MODEL %s_%s\n./_SH/%s_%s_%s.sh\n' % (ident, i, ident, i, args.step))
        f_batch.close()
        subprocess.call(['chmod +rx STEP_%s.sh' % args.step], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

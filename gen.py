import argparse
import os
import shutil
import datetime

### Script argument pares

parser = argparse.ArgumentParser(description='MDGen')
parser.add_argument('-start', action='store_true')
parser.add_argument('-wd', help='Working dir')
parser.add_argument('-top', help='Topology')
parser.add_argument('-i', help='Input for this step')
parser.add_argument('-id', help='ID')
parser.add_argument('-s', help='Step number', type=int)
parser.add_argument('-crd', help='Starting coordinates')
parser.add_argument('-nox', action='store_true',
                    help='Do not write trajectory for this step - usually for minimization or heating')
parser.add_argument('-exe', default='pmemd.MPI', help='Executable name')
parser.add_argument('-ref', action='store_true', help='Ref structure')
parser.add_argument('-ntwprt', type=int, help='Ref structure')
args = parser.parse_args()

### Amber argument parser

amparser = argparse.ArgumentParser(description='MDGen')
amparser.add_argument('-O', action='store_true')
amparser.add_argument('-i')
amparser.add_argument('-o')
amparser.add_argument('-p')
amparser.add_argument('-c')
amparser.add_argument('-r')
amparser.add_argument('-x')
amparser.add_argument('-inf')
amparser.add_argument('-ref')
amparser.add_argument('-np')

if args.start:
    if os.path.isfile(args.top) and os.path.isfile(args.crd) and os.path.isfile(args.i):

        ####### MAKE DIRS AND COPY FILES
        try:
            os.mkdir(args.wd)
        except OSError:
            pass
        try:
            os.mkdir('%s/_1/' % (args.wd))
        except OSError:
            pass
        try:
            os.mkdir('%s/_SH/' % (args.wd))
        except OSError:
            pass
        #### Copy input, topology and starting coordinates
        shutil.copy(args.i, '%s/_1/%s_1.in' % (args.wd, args.id))
        shutil.copy(args.top, '%s/%s_i_s.parm7' % (args.wd, args.id))
        shutil.copy(args.crd, '%s/%s_i_s.rst7' % (args.wd, args.id))
        print ("%s - step 1: Preparing inputs. " % (args.id))
        ####### FILENAMES
        i = '_1/%s_1.in' % (args.id)  # Input
        top = '%s_i_s.parm7' % (args.id)  # Topology
        crd = '%s_i_s.rst7' % (args.id)  # Coordinates
        out = '_1/%s_1.out' % (args.id)  # Output
        inf = '_1/%s_1.info' % (args.id)  # Info
        rst = '_1/%s_1.rst7' % (args.id)  # Restart
        x = '_1/%s_1.nc' % (args.id)  # Trajectory
        f = open('%s/%s.mdg' % (args.wd, args.id), 'w')
        ####### WRITE DATA
        now = datetime.datetime.now()
        date = now.strftime("%Y-%m-%d %H:%M")
        f.write('#MDGen\n')
        f.write(args.id + "\n")
        cmd = '%s -O -i %s -o %s -p %s -c %s -r %s -inf %s' % (args.exe, i, out, top, crd, rst, inf)
        if not args.nox:
            cmd += ' -x %s' % (x)
        f.write('1,%s\n' % (date))
        f.close()
        f = open('%s/_SH/1.sh' % (args.wd), 'w')
        f.write(cmd + "\n")
        f.close()
        print ("OK")
else:
    if os.path.isfile(args.i):
        mdg_files = [f for f in os.listdir(args.wd) if f.endswith('.mdg')]
        f = open('%s/%s' % (args.wd, mdg_files[0]))
        lines = f.readlines()
        ident = args.id
        last_step = int(lines[len(lines) - 1].split(',')[0])
        current_step = last_step + 1
        if args.s:
            current_step = args.s
        print ("%s - step %s: preparing inputs." % (ident, current_step))
        try:
            os.mkdir('%s/_%s/' % (args.wd, current_step))
        except OSError:
            pass
        shutil.copy(args.i, '%s/_%s/%s_%s.in' % (args.wd, current_step, ident, current_step))
        ####### FILENAMES
        amb_args, unk = amparser.parse_known_args(lines[len(lines) - 1].split(',')[1].split()[1:])
        i = '_%s/%s_%s.in' % (current_step, ident, current_step)
        top = '%s_i_s.parm7' % (ident)
        crd = '%s' % (amb_args.r)
        if args.s:
            previous_step = int(current_step) - 1
            crd = '_%s/%s_%s.rst7' % (previous_step, ident, previous_step)
        out = '_%s/%s_%s.out' % (current_step, ident, current_step)
        inf = '_%s/%s_%s.info' % (current_step, ident, current_step)
        rst = '_%s/%s_%s.rst7' % (current_step, ident, current_step)
        x = '_%s/%s_%s.nc' % (current_step, ident, current_step)
        now = datetime.datetime.now()
        date = now.strftime("%Y-%m-%d %H:%M")
        f = open('%s/%s.mdg' % (args.wd, ident), 'a')
        cmd = '%s -O -i %s -o %s -p %s -c %s -r %s -inf %s' % (args.exe, i, out, top, crd, rst, inf)
        if not args.nox:
            cmd += ' -x %s' % (x)
        if args.ref:
            cmd += ' -ref %s' % (crd)
        f.write('%s,%s\n' % (current_step, date))
        f.close()
        if args.ntwprt:
            f = open('%s/_%s/%s_%s.in' % (args.wd, current_step, ident, current_step), 'r')
            lines = f.readlines()
            is_ntwprt = 0
            new_lines = []
            for i in range(0, len(lines)):
                if 'ntwprt' in lines[i]:
                    is_ntwprt += 1
                    cmds = lines[i].split(',')
                    for j in range(0, len(cmds) - 1):
                        if 'ntwprt' in cmds[j]:
                            cmds[j] = 'ntwprt=%s' % (args.ntwprt)
                    new_lines.append(','.join(cmds))

                else:
                    new_lines.append(lines[i])
            if is_ntwprt == 0:
                z = new_lines[3].rstrip()
                z += 'ntwprt=%s,\n' % (args.ntwprt)
                new_lines[3] = z
            f.close()
            f = open('%s/_%s/%s_%s.in' % (args.wd, current_step, ident, current_step), 'w')
            for line in new_lines:
                f.write(line)
            f.close()
        f = open('%s/_SH/%s.sh' % (args.wd, current_step), 'w')
        f.write(cmd + "\n")
        f.close()
        print ("OK")

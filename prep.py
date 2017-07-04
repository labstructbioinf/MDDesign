import argparse
import subprocess
import os
parser = argparse.ArgumentParser(description='MDGen')
parser.add_argument('-pdb')
parser.add_argument('-d', type=float, default=8.0)
parser.add_argument('-ions', action='store_true')
parser.add_argument('-parm')
parser.add_argument('-sparm')
parser.add_argument('-crd')
args = parser.parse_args()
ident = os.path.splitext(os.path.basename(args.pdb))[0]
f = open("%s_leap.in" % (ident), 'w')
inp = (
"""
source leaprc.protein.ff14SB
source leaprc.water.tip3p
mol = loadpdb %s_temp.pdb
solvateoct mol TIP3PBOX %s
""" % (ident, args.d))
if args.ions:
	inp += (
"""addions mol Na+ 0
addions mol Cl- 0""")
inp += (
"""
saveamberparm mol %s %s
quit""" % (args.parm, args.crd))
f.write(inp)
f.close()
subprocess.call(['$AMBERHOME/bin/reduce -Trim %s > %s_temp.pdb' % (args.pdb, ident)], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
subprocess.call(['$AMBERHOME/bin/tleap -f %s_leap.in > %s_leap.log' % (ident, ident)], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
subprocess.call(['rm %s_leap.in %s_leap.log %s_temp.pdb' % (ident, ident, ident)], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

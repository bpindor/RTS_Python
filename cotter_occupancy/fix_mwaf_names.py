import glob
import sys, os


mwaf_files = glob.glob('*.mwaf')

for f in mwaf_files:
    cmd = 'cp %s new_1125953008_%s' % (f, f[-7:])
    os.system(cmd)
    print cmd

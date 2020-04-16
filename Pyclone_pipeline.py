

import subprocess
from glob import glob
import sys


YAMLDIR = '/data/haghirebrahimm2/Pyclone_AA_Lung'
CONFIG_FILE = '/home/haghirebrahimm2/LungCancer/config.yaml'
# TSVDIR = glob(sys.argv[1] + '/*.tsv')



# for TSV in TSVDIR:
# 	PATIENTID = TSV[TSV.rfind('/') + 1: TSV.index('.')]
# 	cmd1 = 'PyClone build_mutations_file --in_file {inp} --out_file {out}'.format(inp = TSV, out = YAMLDIR + '/{}.yaml'.format(PATIENTID))

# 	subprocess.check_output(cmd1, stderr=subprocess.STDOUT, shell = True)


## Set the prior flagk

cmd2 = 'PyClone run_analysis --config_file {conf} --prior  {pri}'.format(conf = CONFIG_FILE, pri = 'total_copy_number')
subprocess.check_output(cmd2, stderr=subprocess.STDOUT, shell = True)
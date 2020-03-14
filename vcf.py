

import sys
import numpy as np
import subprocess
from time import time, sleep
import os
import pandas as pd
from glob import glob

VCFDIR = glob(sys.argv[1] + '/*.FINAL.vcf')
SELECTED_VCFs = []
for vcf in VCFDIR:
	cmd = 'wc -l {}'.format(vcf)
	out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell = True)
	out = out.decode("utf-8").splitlines()
	if int(out[0].split(' ')[0]) < 10000:
		SELECTED_VCFs.append(out[0].split(' ')[1])
	
print(len(SELECTED_VCFs))
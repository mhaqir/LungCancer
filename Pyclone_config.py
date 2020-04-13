

from glob import glob
import sys

WORKING_DIR = '/data/haghirebrahimm2/Pyclone_AA_Lung'
TRACE_DIR = '/data/haghirebrahimm2/Pyclone_AA_Lung/Pylone_output'
DENSITY = 'pyclone_beta_binomial'
NUM_ITERS = '10000'
BASE_MEASURE_PARAM_ALPHA = '1'
BASE_MEASURE_PARAM_BETA = '1'
CONCENTRATION_VALUE = '1.0'
PRIOR_SHAPE = '1.0'
PRIOR_RATE = '0.001'
SAMPLE_ID = 'AfAm'

TSVDIR = glob(sys.argv[1] + '/*.tsv')



with open('config.txt', 'w') as f:
	f.writelines('working_dir: {}\n'.format(WORKING_DIR))
	f.writelines('trace_dir: {}\n'.format(TRACE_DIR))
	f.writelines('density: {}\n'.format(NUM_ITERS))
	f.writelines('base_measure_params:\n')
	f.writelines('	alpha: {}\n'.format(BASE_MEASURE_PARAM_ALPHA))
	f.writelines('	beta: {}\n'.format(BASE_MEASURE_PARAM_BETA))
	f.writelines('concentration:\n')
	f.writelines('	value: {}\n'.format(CONCENTRATION_VALUE))
	f.writelines('prior:\n')
	f.writelines('	shape: {}\n'.format(PRIOR_SHAPE))
	f.writelines('	rate: {}\n'.format(PRIOR_RATE))
	f.writelines('samples:\n')
	f.writelines('	{}:\n'.format(SAMPLE_ID))
	for i in range(len(TSVDIR) - 1):
		f.writelines('		mutations_file: {}\n'.format(TSVDIR[i]))
	f.writelines('		mutations_file: {}'.format(TSVDIR[len(TSVDIR) - 1]))


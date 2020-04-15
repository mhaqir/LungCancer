

from glob import glob
import sys

WORKING_DIR = '/data/haghirebrahimm2/Pyclone_AA_Lung'
TRACE_DIR = 'Pyclone_output/trace'
DENSITY = 'pyclone_beta_binomial'
NUM_ITERS = '10000'
BASE_MEASURE_PARAM_ALPHA = '1'
BASE_MEASURE_PARAM_BETA = '1'
CONCENTRATION_VALUE = '1.0'
PRIOR_SHAPE = '1.0'
PRIOR_RATE = '0.001'
SAMPLE_ID = 'AfAm'
TUMOUR_CONTENT_VALUE = '1.0'
ERROR_RATE = '0.001'
BETA_BINOMIAL_PRECISION_PARAMS = {'value': '1000', 'prior_shape': '1.0',\
 'prior_rate': '0.0001', 'proposal_precision': '0.01'}

YAMLDIR = glob(sys.argv[1] + '/*.yaml')



with open('config.yaml', 'w') as f:
	f.writelines('num_iters: {}\n\n'.format(NUM_ITERS))

	f.writelines('base_measure_params:\n')
	f.writelines('  alpha: {}\n'.format(BASE_MEASURE_PARAM_ALPHA))
	f.writelines('  beta: {}\n\n'.format(BASE_MEASURE_PARAM_BETA))

	f.writelines('concentration:\n')
	f.writelines('  value: {}\n\n'.format(CONCENTRATION_VALUE))
	f.writelines('  prior:\n')
	f.writelines('    shape: {}\n'.format(PRIOR_SHAPE))
	f.writelines('    rate: {}\n\n'.format(PRIOR_RATE))

	f.writelines('density: {}\n\n'.format(DENSITY))

	f.writelines('beta_binomial_precision_params:\n')
	f.writelines('  value: {}\n\n'.format(BETA_BINOMIAL_PRECISION_PARAMS['value']))
	f.writelines('  prior:\n')
	f.writelines('    shape: {}\n'.format(BETA_BINOMIAL_PRECISION_PARAMS['prior_shape']))
	f.writelines('    rate: {}\n\n'.format(BETA_BINOMIAL_PRECISION_PARAMS['prior_rate']))
	f.writelines('  proposal:\n')
	f.writelines('    precision: {}\n\n'.format(BETA_BINOMIAL_PRECISION_PARAMS['proposal_precision']))

	f.writelines('working_dir: {}\n\n'.format(WORKING_DIR))

	f.writelines('trace_dir: {}\n\n'.format(TRACE_DIR))

	f.writelines('samples:\n')
	f.writelines('  {}:\n'.format(SAMPLE_ID))
	for i in range(len(YAMLDIR) - 1):
		YAMLFILE = YAMLDIR[i][YAMLDIR[i].rfind('/') + 1:]
		f.writelines('    mutations_file: {}\n'.format(YAMLFILE))
	YAMLFILE = YAMLDIR[len(YAMLDIR) - 1][YAMLDIR[len(YAMLDIR) - 1].rfind('/') + 1:]
	f.writelines('    mutations_file: {}\n\n'.format(YAMLFILE))
	f.writelines('    tumour_content:\n')
	f.writelines('      value: {}\n\n'.format(TUMOUR_CONTENT_VALUE))
	f.writelines('    error_rate: {}'.format(ERROR_RATE))

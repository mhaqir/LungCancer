#!/bin/python3


from glob import glob
import sys
import numpy as np
import pandas as pd
from time import time

from sklearn.metrics.pairwise import cosine_distances
from sklearn.neighbors import DistanceMetric


from skbio import DistanceMatrix
from skbio.tree import nj

VCFDIR = glob(sys.argv[1] + '/*.FINAL.vcf')
CHROMS = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', \
'12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']

def records(VCF):
	with open(VCF, 'r') as f:
		lines = f.readlines()

	lines = [item for item in lines if item[0] != '#']    # Removing meta-information
	SNV_records = [item.rstrip().split('\t') for item in lines[1:]]    # Extracting SNV records

	muts = []
	# c = 0
	# c1 = 0
	for record in SNV_records:
		if record[0] in ['X', 'Y']:
			# print('x')
			continue
		if record[6] != 'PASS':
			# print('xx')
			continue

		INFO = record[7].split(';')
		if INFO[-1][0:2] != 'VT':
			# c += 1
			# print('xxx')
			continue

		# FORMAT = record[8]
		# if FORMAT != 'GT:AD:AF':
		# 	continue
		FORMAT = record[8]
		if FORMAT != 'GT:AD:BQ:DP:FA:SS':
			# print('xxxx')
			continue

		if record[0] not in CHROMS:
			# print('xxxxx')
			continue

		# r7.append(record[7])
		muts.append('chr' + record[0] + ':' + record[1])
		# c1 += 1
	return muts


mutations = {}
PIDS = []
for VCF in VCFDIR:
	PATIENTID = VCF[VCF.index('+') + 1:VCF.index('.')]
	PIDS.append(PATIENTID)
	mutations[PATIENTID] = records(VCF)


mutations_list = [v for k,v in mutations.items()]
mutations_list = list(set([item[i] for item in mutations_list for i in range(len(item))]))

samples_mutations = pd.DataFrame(np.zeros((len(VCFDIR), len(mutations_list)), dtype = int), columns = mutations_list, \
									index = PIDS)




store = pd.HDFStore('samples_mutations.h5')
store['samples_mutations'] = samples_mutations

# samples_mutations_load = store['samples_mutations']




# dm_cosine = DistanceMatrix(cosine_distances(samples_mutations_load.values), PIDS)
# # print(dm)

# tree_cosine = nj(dm_cosine)

# print('Tree based on Cosine distance:')
# print(tree_cosine.ascii_art())



# dist = DistanceMetric.get_metric('hamming')
# dm_hamming = DistanceMatrix(dist.pairwise(samples_mutations_load.values), PIDS)

# tree_hamming = nj(dm_hamming)
# print('Tree based on Hamming distance:')
# print(tree_hamming.ascii_art())



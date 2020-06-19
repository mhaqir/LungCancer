

from glob import glob
import subprocess
from os import path
import pandas as pd
import numpy as np
from collections import Counter


VCFDIR_GDC = '/data/haghirebrahimm2/GDC_LungCancer'
VCFDIR_BMR = '/data/BMR_Genomics/data_for_gardner/exomeseq_activeDev/mutect_out'

CLINICALDIR_GDC = '/home/haghirebrahimm2/LungCancer/GDC_clinical'

def records_GDC(VCF, chrs):
	with open(VCF, 'r') as f:
		lines = f.readlines()
	# print(len(lines))
	lines = [item for item in lines if item[0] != '#']    # Removing meta-information
	SNV_records = [item.rstrip().split('\t') for item in lines[1:]]    # Extracting SNV records

	muts = []
	# c = 0
	# c1 = 0
	for record in SNV_records:
		if record[0] in ['chrX', 'chrY']:
			# print('x')
			continue
		if record[6] != 'PASS':
			# print('xx')
			continue

		# INFO = record[7].split(';')
		# if INFO[-1][0:2] != 'VT':
			# c += 1
			# print('xxx')
			# continue

		# FORMAT = record[8]
		# if FORMAT != 'GT:AD:AF':
		# 	continue
		# FORMAT = record[8]
		# if FORMAT != 'GT:AD:BQ:DP:FA:SS':
			# print('xxxx')
			# continue

		if record[0] not in chrs:
			# print('xxxxx')
			continue

		# r7.append(record[7])
		muts.append(record[0] + ':' + record[1])
		# c1 += 1
	return muts


def records_BMR(VCF, chrs):
	with open(VCF, 'r') as f:
		lines = f.readlines()
	# print(len(lines))
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

		# INFO = record[7].split(';')
		# if INFO[-1][0:2] != 'VT':
			# c += 1
			# print('xxx')
			# continue

		# FORMAT = record[8]
		# if FORMAT != 'GT:AD:AF':
		# 	continue
		# FORMAT = record[8]
		# if FORMAT != 'GT:AD:BQ:DP:FA:SS':
			# print('xxxx')
			# continue

		if record[0] not in chrs:
			# print('xxxxx')
			continue

		# r7.append(record[7])
		muts.append('chr' + record[0] + ':' + record[1])
		# c1 += 1
	return muts




###############  GDC  ###################
foldernames_GDC = glob(VCFDIR_GDC + '/*')

downloaded_samples = [f[f.rfind('/') + 1:] for f in foldernames_GDC]

clinical = pd.read_csv(CLINICALDIR_GDC + '/clinical.tsv', sep = '\t')
sample_sheet = pd.read_csv(CLINICALDIR_GDC + '/gdc_sample_sheet.tsv', sep = '\t')

sample_sheet_case_id = [item.split(', ')[0] for item in list(sample_sheet['Case ID'].values)]
sample_sheet_file_id = list(sample_sheet['File ID'].values)
assert len(sample_sheet_file_id) == len(sample_sheet_case_id)
sample_sheet_dict = {sample_sheet_file_id[i]: sample_sheet_case_id[i] for i in range(len(sample_sheet_file_id))}


c = Counter(sample_sheet_case_id)
downloaded_samples_unique_case_id = {s: sample_sheet_dict[s] for s in downloaded_samples if c[sample_sheet_dict[s]] == 1}


samples_table = []
for k, v  in downloaded_samples_unique_case_id.items():
	row = clinical.iloc[np.where(clinical['case_submitter_id'].values == v)[0], ]
	if (np.shape(row.values)[0] > 0) and (row['race'].values[0] == 'black or african american') and (row['project_id'].values[0] in ['TCGA-LUSC', 'TCGA-LUAD']):
		samples_table.append([k, v, row['race'].values[0], row['gender'].values[0], row['project_id'].values[0]])

samples_table = np.asarray(samples_table)

CHROMS = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', \
'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']


mutations_GDC = {}
PIDS_GDC = []
for i in range(np.shape(samples_table)[0]):
	PATIENTID = samples_table[i, 1]
	PIDS_GDC.append(PATIENTID)
	mutations_GDC[PATIENTID] = records_GDC(VCFDIR_GDC + '/' + samples_table[i, 0] + '/' + samples_table[i, 0] + '.vcf', CHROMS)
############### BMR #################

filenames_BMR = glob(VCFDIR_BMR + '/*.FINAL.vcf')

CHROMS = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', \
'12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']

mutations_BMR = {}
PIDS_BMR = []
for VCF in filenames_BMR:
	PATIENTID = VCF[VCF.index('+') + 1:VCF.index('.')]
	PIDS_BMR.append(PATIENTID)
	mutations_BMR[PATIENTID] = records_BMR(VCF, CHROMS)
##################################

mutations_list_GDC = [v for k,v in mutations_GDC.items()]
mutations_list_BMR = [v for k,v in mutations_BMR.items()]

mutations_list = mutations_list_GDC + mutations_list_BMR

# print(len(mutations_list_BMR))
# print(len(mutations_list_GDC))
# print(len(mutations_list))

mutations_list = list(set([item[i] for item in mutations_list for i in range(len(item))]))
PIDS = PIDS_GDC + PIDS_BMR
mutations = {}
mutations.update(mutations_GDC)
mutations.update(mutations_BMR)

# print(len(mutations))

samples_mutations = pd.DataFrame(np.zeros((len(PIDS), len(mutations_list)), dtype = int), columns = mutations_list, \
									index = PIDS)

for k, v in mutations.items():
	for mut in v:
		samples_mutations.loc[k, mut] = 1



store = pd.HDFStore('GDCbaa_BMR_samples_mutations.h5')
store['GDCbaa_BMR_samples_mutations'] = samples_mutations


from glob import glob
import subprocess
from os import path
import pandas as pd
import numpy as np
from collections import Counter

from GTF_preproc import readGTF, gene

VCFDIR = '/data/haghirebrahimm2/GDC_LungCancer'
CLINICALDIR = '/home/haghirebrahimm2/LungCancer/GDC_clinical'
GTFDIR = '/home/haghirebrahimm2/refgenome/hg19refGene.gtf.sorted'
CANCER_CENSUS_GENE_DIR = '/home/haghirebrahimm2/LungCancer/cancer_gene_census.csv'

foldernames = glob(VCFDIR + '/*')

# for f in foldernames:
# 	if path.exists(f + '/' + f[f.rfind('/') + 1:] + '.vcf'):
# 		continue
# 	else:
# 		vcf = f + '/' + f[f.rfind('/') + 1:] + '.vcf.gz'
# 		subprocess.call(['gzip', '-d', vcf])


Cancer_census_genes = pd.read_csv(CANCER_CENSUS_GENE_DIR, sep =',')['Gene Symbol'].values
downloaded_samples = [f[f.rfind('/') + 1:] for f in foldernames]

clinical = pd.read_csv(CLINICALDIR + '/clinical.tsv', sep = '\t')
sample_sheet = pd.read_csv(CLINICALDIR + '/gdc_sample_sheet.tsv', sep = '\t')

sample_sheet_case_id = [item.split(', ')[0] for item in list(sample_sheet['Case ID'].values)]
sample_sheet_file_id = list(sample_sheet['File ID'].values)
assert len(sample_sheet_file_id) == len(sample_sheet_case_id)
sample_sheet_dict = {sample_sheet_file_id[i]: sample_sheet_case_id[i] for i in range(len(sample_sheet_file_id))}


c = Counter(sample_sheet_case_id)
downloaded_samples_unique_case_id = {s: sample_sheet_dict[s] for s in downloaded_samples if c[sample_sheet_dict[s]] == 1}


samples_table = []
for k, v  in downloaded_samples_unique_case_id.items():
	row = clinical.iloc[np.where(clinical['case_submitter_id'].values == v)[0], ] #
	if (np.shape(row.values)[0] > 0) and (row['race'].values[0] in ['black or african american', 'white']) and (row['project_id'].values[0] in ['TCGA-LUSC', 'TCGA-LUAD']):
		samples_table.append([k, v, row['race'].values[0], row['gender'].values[0], row['project_id'].values[0]])


with open('GDC_samples_white_BAA.txt', 'w') as f:
	f.writelines('	'.join(['file_id', 'case_id', 'race', 'gender', 'project_id']) + '\n')
	for i in range(len(samples_table) - 1):
		f.writelines('	'.join(samples_table[i]) + '\n')
	f.writelines('	'.join(samples_table[-1]))

samples_table = np.asarray(samples_table)

CHROMS = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', \
'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

def records(VCF, GTF_dict):
	with open(VCF, 'r') as f:
		lines = f.readlines()
	# print(len(lines))
	lines = [item for item in lines if item[0] != '#']    # Removing meta-information
	SNV_records = [item.rstrip().split('\t') for item in lines[1:]]    # Extracting SNV records

	muts = []
	# c = 0
	# c1 = 0
	for record in SNV_records:
		# if record[0] in ['chrX', 'chrY']:
			# continue
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

		if record[0] not in CHROMS:
			# print('xxxxx')
			continue

		# r7.append(record[7])
		gene_id = gene(record[0][3:], int(record[1]), GTF_dict)
		if (gene_id != '') and (gene_id in Cancer_census_genes):
			muts.append(record[0] + ':' + record[1])
		# c1 += 1
	return muts

GTF_dict = readGTF(GTFDIR)
mutations = {}
PIDS = []
for i in range(np.shape(samples_table)[0]):
	PATIENTID = samples_table[i, 1]
	muts = records(VCFDIR + '/' + samples_table[i, 0] + '/' + samples_table[i, 0] + '.vcf', GTF_dict)
	if len(muts) > 0:
		PIDS.append(PATIENTID)
		mutations[PATIENTID] = muts


mutations_list = [v for k,v in mutations.items()]
mutations_list = list(set([item[i] for item in mutations_list for i in range(len(item))]))

samples_mutations = pd.DataFrame(np.zeros((len(PIDS), len(mutations_list)), dtype = int), columns = mutations_list, \
									index = PIDS)

for k, v in mutations.items():
	for mut in v:
		samples_mutations.loc[k, mut] = 1


store = pd.HDFStore('GDCw_baa_samples_mutations_cc.h5')
store['GDCw_baa_samples_mutations_cc'] = samples_mutations
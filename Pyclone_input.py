


import sys
from glob import glob
import pandas as pd

OUTPUTDIR = '/data/haghirebrahimm2/Pyclone_AA_Lung'
VCFDIR = glob(sys.argv[1] + '/*.FINAL.vcf')


clinical = pd.read_excel('Clinical.xlsx', 'Clinical Data')
samples_excel = list(clinical['Tumor_Sample_Barcode'].values)
samples_gender = list(clinical['Gender'].values)


CHROMS = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', \
'12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']


def records(VCF, PATIENTID):
	with open(VCF, 'r') as f:
		lines = f.readlines()

	lines = [item for item in lines if item[0] != '#']    # Removing meta-information
	SNV_records = [item.rstrip().split('\t') for item in lines[1:]]    # Extracting SNV records

	Pyclone_records = []
	gender =  samples_gender[samples_excel.index(PATIENTID)] ##'unknown'  ##sys.argv[2]
	for record in SNV_records:

		if record[6] != 'PASS':
			continue

		INFO = record[7].split(';')
		if INFO[1][0:2] != 'VT':
			continue

		# FORMAT = record[8]
		# if FORMAT != 'GT:AD:AF':
		# 	continue
		FORMAT = record[8]
		if FORMAT != 'GT:AD:BQ:DP:FA:SS':
			continue

		if record[0] not in CHROMS:
			continue

		normal = record[9]
		tumor = record[10].split(':')
		tumor_AD_values = tumor[1].split(',')
		tumor_REF_cov = tumor_AD_values[0]
		tumor_ALT_cov = tumor_AD_values[1]

		if record[0] in ['X', 'Y']:
			normal_cn = 1
		else:
			normal_cn = 2

		minor_cn = 0
		major_cn = 2

		Pyclone_records.append([record[0] + '_' + record[1], tumor_REF_cov, tumor_ALT_cov, str(normal_cn), str(minor_cn), str(major_cn)])
	return Pyclone_records



for VCF in VCFDIR:
	PATIENTID = VCF[VCF.index('+') + 1:VCF.index('.')]
	Pyclone_input = records(VCF, PATIENTID)

	Pyclone_col_names = ['mutation_id', 'ref_counts', 'var_counts','normal_cn', 'minor_cn', 'major_cn']
	with open(OUTPUTDIR + '/' + '{}.tsv'.format(PATIENTID), 'w') as f:
		f.writelines('	'.join(Pyclone_col_names) + '\n')
		for i in range(len(Pyclone_input) - 1):
			f.writelines('	'.join(Pyclone_input[i]) + '\n')
		f.writelines('	'.join(Pyclone_input[len(Pyclone_input)-1]))




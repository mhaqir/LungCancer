
import sys


with open(sys.argv[1], 'r') as f:
	lines = f.readlines()

lines = [item for item in lines if item[0] != '#']    # Removing meta-information
# col_names = lines[0].rstrip().split('\t')
# col_names[0] = col_names[0][1:]    # Extracting feild names
SNV_records = [item.rstrip().split('\t') for item in lines[1:]]    # Extracting SNV records



CTPsingle_record = []
multiplier = 2
gender = sys.argv[2]

for record in SNV_records:
	# if record[0] in ['X', 'Y']:
	# 	continue
	if record[6] != 'PASS':
		continue

	INFO = record[7].split(';')
	if INFO[4] != 'SOMATIC':
		continue

	FORMAT = record[8]
	if FORMAT != 'GT:AD:AF':
		continue
	normal = record[9]
	tumor = record[10].split(':')
	tumor_AD_values = tumor[1].split(',')
	tumor_REF_cov = tumor_AD_values[0]
	tumor_ALT_cov = tumor_AD_values[1]

	CTPsingle_record.append([record[0], record[1], record[4], record[3], tumor_ALT_cov, tumor_REF_cov, str(multiplier), gender])


CTPsingle_col_names = ['Chromosome', 'Position', 'Mutatnt', 'Refrence', 'Mcount', 'Rcount', 'Multiplier', 'Gender']
filename = sys.argv[1][sys.argv[1].index('+') + 1:sys.argv[1].index('.')]

with open(filename + '.frq', 'w') as f:
	f.writelines('\t'.join(CTPsingle_col_names) + '\n')
	for record in CTPsingle_record:
		f.writelines('\t'.join(record) + '\n')





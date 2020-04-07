
import sys
import numpy as np
import subprocess
from time import time, sleep
import os
import pandas as pd
from glob import glob


from GTF_preproc import readGTF, gene



# (chrom, position, VAF, cluster, alteration_type, gene_id)
CANCER_CENSUS_GENE_DIR = '/home/haghirebrahimm2/LungCancer/cancer_gene_census.csv'
GTFDIR = '/home/haghirebrahimm2/refgenome/hg19refGene.gtf.sorted'
OUTPUTDIR = '/data/haghirebrahimm2/CTPsingle_AA_Lung'

Cancer_census_genes = pd.read_csv(CANCER_CENSUS_GENE_DIR, sep =',')['Gene Symbol'].values
VCFDIR = glob(sys.argv[1] + '/*.FINAL.vcf')
SELECTED_VCFs = []
for vcf in VCFDIR:
	cmd = 'wc -l {}'.format(vcf)
	out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell = True)
	out = out.decode("utf-8").splitlines()
	if int(out[0].split(' ')[0]) < 10000:
		SELECTED_VCFs.append(out[0].split(' ')[1])

clinical = pd.read_excel('Clinical.xlsx', 'Clinical Data')
samples_excel = list(clinical['Tumor_Sample_Barcode'].values)
samples_gender = list(clinical['Gender'].values)

CHROMS = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', \
'12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']

def records(VCF):
	with open(VCF, 'r') as f:
		lines = f.readlines()

	lines = [item for item in lines if item[0] != '#']    # Removing meta-information
	SNV_records = [item.rstrip().split('\t') for item in lines[1:]]    # Extracting SNV records

	CONETT_records = []
	CTPsingle_records = []
	multiplier = 2
	gender =  samples_gender[samples_excel.index(PATIENTID)] ##'unknown'  ##sys.argv[2]
	for record in SNV_records:
		if record[0] in ['X', 'Y']:
			continue
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

		# CTPsingle_record.append([record[0], record[1], record[4], record[3], tumor_ALT_cov, tumor_REF_cov,\
		#  (2 * tumor_ALT_cov)/tumor_REF_cov, INFO[1][3:]])

		CONETT_records.append([record[0], record[1], (2 * int(tumor_ALT_cov))/(int(tumor_REF_cov) + int(tumor_ALT_cov)), INFO[1][3:]])
		CTPsingle_records.append([record[0], record[1], record[4], record[3], tumor_ALT_cov, tumor_REF_cov, str(multiplier), gender])
	return CTPsingle_records, CONETT_records

# VCF = sys.argv[1] + '/combined_SC284874+combined_SC284875.FINAL.vcf' 
# VCF = VCFDIR[0]

# Creating dict from GTF file and extracting gene names
GTF_dict = readGTF(GTFDIR)
CONETT_input = []

for VCF in SELECTED_VCFs:
	PATIENTID = VCF[VCF.index('+') + 1:VCF.index('.')]
	if not os.path.exists(OUTPUTDIR + '/' + PATIENTID):
		os.mkdir(OUTPUTDIR + '/' + PATIENTID)

		CTPsingle_records, CONETT_records = records(VCF)
		# CTPsingle input
		CTPsingle_col_names = ['Chromosome', 'Position', 'Mutant', 'Refrence', 'Mcount', 'Rcount', 'Multiplier', 'Gender']
		# filename = sys.argv[1][sys.argv[1].index('+') + 1:sys.argv[1].index('.')]

		with open(OUTPUTDIR + '/' + PATIENTID + '/' + PATIENTID + '.frq', 'w') as f:
			f.writelines(' '.join(CTPsingle_col_names) + '\n')
			for i in range(len(CTPsingle_records) - 1):
				f.writelines(' '.join(CTPsingle_records[i]) + '\n')
			f.writelines(' '.join(CTPsingle_records[len(CTPsingle_records)-1]))

		# Run CTPsingle
		# cmd = ['/home/haghirebrahimm2/Tools/CTPsingle/CTPsingle.R','-f', OUTPUTDIR + '/' + PATIENTID + '/' + PATIENTID + '.frq',\
		# '-o', OUTPUTDIR + '/' + PATIENTID + '/' + PATIENTID, '-m', '/home/haghirebrahimm2/Tools/CTPsingle/GammaAdjMatrices']
		# print('1')
		# subprocess.run(cmd)
		# print('2')

		cmd = '/home/haghirebrahimm2/Tools/CTPsingle/CTPsingle.R ' + '-f ' + OUTPUTDIR + '/' + PATIENTID + '/' + PATIENTID + '.frq ' +\
		'-o ' + OUTPUTDIR + '/' + PATIENTID + '/' + PATIENTID + ' -m ' + '/home/haghirebrahimm2/Tools/CTPsingle/GammaAdjMatrices'
		# os.system(cmd)
		subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell = True)

		## CTP_single output
		with open(OUTPUTDIR + '/' + PATIENTID + '/' + PATIENTID + '_cluster_assignments.txt', 'r') as f:
			lines = f.readlines()

		CTP_output = [item.rstrip().split('\t') for item in lines[1:]]
		CTP_output = np.asarray(CTP_output)
		CTP_Clusters = CTP_output[:,4]

		CONETT_records = np.asarray(CONETT_records)
		CTP_Clusters = np.expand_dims(np.asarray(CTP_Clusters), axis = 1)

		# CONETT_records = CONETT_records[CONETT_records[:, 0:2] == CTP_output[:, 0:2]]

		# Appending clusters from CTPsingle
		assert np.shape(CONETT_records)[0] == np.shape(CTP_Clusters)[0]
		CONETT_records = np.concatenate((CONETT_records, CTP_Clusters), axis = 1)

		genes = []
		idx = []
		# st = time()
		# Extracting gene ids and apply cancer-genes filter
		for i in range(np.shape(CONETT_records)[0]):
			gene_id = gene(CONETT_records[i,0], int(CONETT_records[i,1]), GTF_dict)
			if (gene_id != '') and (gene_id in Cancer_census_genes):
				genes.append(gene_id)
				idx.append(i)
		# print(time() - st)

		genes = np.expand_dims(np.asarray(genes), axis = 1)
		CONETT_records = CONETT_records[idx, :]

		assert np.shape(CONETT_records)[0] == np.shape(genes)[0]
		CONETT_records = np.concatenate((CONETT_records, genes), axis = 1)

		# Cluster labels and their mean VAF
		Cluster_labels = list(set(list(CONETT_records[:,4])))
		# print(Cluster_labels)

		# Clusters mean VAF
		# Clutser_VAFS = [np.mean(CONETT_records[CONETT_records[:, 4] == l][:, 2].astype(float)) for l in Cluster_labels]

		# In case there are more than one mutation on a gene, we get the one with highest VAF 
		CONETT_selected_records = []
		for l in Cluster_labels:
			group = CONETT_records[CONETT_records[:, 4] == l, :] # Group by class
			genes_in_group = list(set(list(group[:, 5])))  # Unique genes in a class
			for gene_name in genes_in_group:
				group_g = group[group[:, 5] == gene_name, :]   # Group by genes in a class
				selected_record = group_g[np.argmax(group_g[:, 2].astype(float)), :]   # Select the record with highest VAF

				CONETT_selected_records.append(list(selected_record))

		# print(len(CONETT_selected_records))
		# print(CONETT_selected_records[0:10])


		CONETT_selected_records = np.asarray(CONETT_selected_records)
		# print(CONETT_selected_records)
		# print(np.shape(CONETT_selected_records))
		Cluster_labels_selected = list(set(list(CONETT_selected_records[:,4])))
		Cluster_VAFS = [np.mean(CONETT_selected_records[CONETT_selected_records[:, 4] == l][:, 2].astype(float)) for l in Cluster_labels_selected]

		# print('Cluster_labels_selected', Cluster_labels_selected)
		# print('Cluster_VAFS', Cluster_VAFS)

		Cluster_pairs = [(c1, c2) for c1 in Cluster_labels_selected for c2 in Cluster_labels_selected \
		 if Cluster_VAFS[Cluster_labels_selected.index(c1)] > Cluster_VAFS[Cluster_labels_selected.index(c2)]]

		# print('Clutser_VAFS', Cluster_pairs)
		
		for pair in Cluster_pairs:
			group1 = CONETT_selected_records[CONETT_selected_records[:, 4] == pair[0]]
			# print('group1 {}'.format(pair), np.shape(group1))
			group2 = CONETT_selected_records[CONETT_selected_records[:, 4] == pair[1]]
			# print('group2 {}'.format(pair), np.shape(group2))
			CONETT_input = CONETT_input + [[PATIENTID, group1[i, 5], group1[i, 3], group1[i, 4], group2[j, 5], group2[j, 3], group2[j, 4]] \
			                 for i in range(np.shape(group1)[0]) for j in range(np.shape(group2)[0])]# \
			                 # if (group1[i, 5] != '') and (group2[j, 5] != '') and \
			                 # (group1[i, 5] in Cancer_census_genes) and (group2[j, 5] in Cancer_census_genes) ]
			# print('CONETT_input', len(CONETT_input))
		# sleep(5)

# print(len(CONETT_input))
# print(CONETT_input[0:10])
CONETT_col_names = ['patient_id', 'gene_name_src', 'alteration_type_src', 'src_clus','gene_name_dest', 'alteration_type_dest', 'des_clus']
with open(OUTPUTDIR + '/' + 'CONETT.txt', 'w') as f:
	f.writelines(' '.join(CONETT_col_names) + '\n')
	for i in range(len(CONETT_input) - 1):
		f.writelines(' '.join(CONETT_input[i]) + '\n')
	f.writelines(' '.join(CONETT_input[len(CONETT_input)-1]))




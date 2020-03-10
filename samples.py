

from glob import glob
# import sys
# from pandas import read_excel
import numpy as np
import pandas as pd


FASTQDIR = '/data/BMR_Genomics/data_for_gardner/exomeseq_activeDev'
MUTECT2DIR = '/data/BMR_Genomics/data_for_gardner/exomeseq_activeDev/mutect2_out' 


Samples_from_fastq_dir = [name[name.rfind('/') + 1:name.index('.')] for name in glob(FASTQDIR + '/*.R1.fastq.gz')]
Samples_from_mutect2_dir = [name[name.rfind('/') + 1:name.index('.')] for name in glob(MUTECT2DIR + '/*.FINALmutect2.vcf')]

clinical = pd.read_excel('Clinical.xlsx', 'Clinical Data')
samples_excel = clinical['Tumor_Sample_Barcode'].values

normal_ids_from_mutect2_dir = [name[0:name.index('+')] for name in Samples_from_mutect2_dir]
tumor_ids_from_mutect2_dir = [name[name.index('+') + 1:] for name in Samples_from_mutect2_dir]



idx = 1
samples = {}
for sample_id in Samples_from_fastq_dir:
	if sample_id not in normal_ids_from_mutect2_dir + tumor_ids_from_mutect2_dir:
		e = [sample_id, sample_id, None]
		if sample_id in samples_excel:
			e.append(1)
		else:
			e.append(None)
		samples[idx] = e
		idx += 1
		continue
	if sample_id in normal_ids_from_mutect2_dir:
		e = [sample_id, tumor_ids_from_mutect2_dir[normal_ids_from_mutect2_dir.index(sample_id)], 1]
		if tumor_ids_from_mutect2_dir[normal_ids_from_mutect2_dir.index(sample_id)] in samples_excel:
			e.append(1)
		else:
			e.append(None)
		samples[idx] = e
		idx += 1
		continue


values = [v for k, v in samples.items()]
keys = [k for k, v in samples.items()]
print(values)
print(keys)
df = pd.DataFrame(values, index = keys, columns = ['fastq_normal', 'fastq_tumor', 'mutect2', 'clinical'])
print(df)
df.to_csv('samples.csv', sep = ',')
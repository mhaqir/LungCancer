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


clinical = pd.read_excel('Clinical.xlsx', 'Clinical Data')
# samples_excel = list(clinical['Tumor_Sample_Barcode'].values)
# samples_gender = list(clinical['Gender'].values)
# samples_histology = list(clinical['Histology'].values)

samples_id_gen_hist = clinical[['Tumor_Sample_Barcode', 'Gender', 'Histology']].values
# print(samples_id_gen_hist)

store = pd.HDFStore('samples_mutations.h5')
samples_mutations_load = store['samples_mutations']

# print(samples_gender)
# print(samples_histology)

# dm_cosine = DistanceMatrix(cosine_distances(samples_mutations_load.values), PIDS)
# tree_cosine = nj(dm_cosine)

# print('Tree based on Cosine distance:')
# print(tree_cosine.ascii_art())

# dist = DistanceMetric.get_metric('hamming')
# dm_hamming = DistanceMatrix(dist.pairwise(samples_mutations_load.values), PIDS)

# tree_hamming = nj(dm_hamming)
# print('Tree based on Hamming distance:')
# print(tree_hamming.ascii_art())

##
# samples_mutations = pd.DataFrame(columns = samples_mutations_load.columns)
# st = time()
# for sample in samples_id_gen_hist[:, 0]:
# 	samples_mutations.loc[sample] = np.squeeze(samples_mutations_load.loc[[sample], :].values)
# print(time() - st)
##


# print((samples_id_gen_hist[:, 0] == np.asarray(samples_mutations.index)))
female_samples = samples_id_gen_hist[(samples_id_gen_hist[:, 1] == 'Female'), 0]
male_samples = samples_id_gen_hist[(samples_id_gen_hist[:, 1] == 'Male'), 0]
# samples_mutations_male = samples_mutations_load[pd.Series(samples_id_gen_hist[:, 1] == 'Male', name = 'bools').values]
# samples_mutations_female = samples_mutations_load[pd.Series(samples_id_gen_hist[:, 1] == 'Female', name = 'bools').values]
samples_mutations_female = pd.DataFrame(columns = samples_mutations_load.columns)
samples_mutations_male = pd.DataFrame(columns = samples_mutations_load.columns)

for sample in male_samples:
	samples_mutations_male.loc[sample] = np.squeeze(samples_mutations_load.loc[[sample], :].values)

for sample in female_samples:
	samples_mutations_female.loc[sample] = np.squeeze(samples_mutations_load.loc[[sample], :].values)


# print(samples_mutations_male.index)
# print(samples_mutations_female.index)

male_cosine = DistanceMatrix(cosine_distances(samples_mutations_male.values), samples_mutations_male.index)
male_tree_cosine = nj(male_cosine)

print('Tree for males based on Cosine distance:')
print(male_tree_cosine.ascii_art())
print('\n\n')


female_cosine = DistanceMatrix(cosine_distances(samples_mutations_female.values), samples_mutations_female.index)
female_tree_cosine = nj(female_cosine)

print('Tree for females based on Cosine distance:')
print(female_tree_cosine.ascii_art())
print('\n\n')


female_LUAD_samples = samples_id_gen_hist[(samples_id_gen_hist[:, 1] == 'Female') & (samples_id_gen_hist[:, 2] == 'LUAD'), 0]
male_LUAD_samples = samples_id_gen_hist[(samples_id_gen_hist[:, 1] == 'Male') & (samples_id_gen_hist[:, 2] == 'LUAD'), 0]
female_LUSC_samples = samples_id_gen_hist[(samples_id_gen_hist[:, 1] == 'Female') & (samples_id_gen_hist[:, 2] == 'LUSC'), 0]
male_LUSC_samples = samples_id_gen_hist[(samples_id_gen_hist[:, 1] == 'Male') & (samples_id_gen_hist[:, 2] == 'LUSC'), 0]


samples_mutations_female_LUAD = pd.DataFrame(columns = samples_mutations_load.columns)
samples_mutations_male_LUAD = pd.DataFrame(columns = samples_mutations_load.columns)
samples_mutations_female_LUSC = pd.DataFrame(columns = samples_mutations_load.columns)
samples_mutations_male_LUSC = pd.DataFrame(columns = samples_mutations_load.columns)


for sample in male_LUAD_samples:
	samples_mutations_male_LUAD.loc[sample] = np.squeeze(samples_mutations_load.loc[[sample], :].values)

for sample in female_LUAD_samples:
	samples_mutations_female_LUAD.loc[sample] = np.squeeze(samples_mutations_load.loc[[sample], :].values)

for sample in male_LUSC_samples:
	samples_mutations_male_LUSC.loc[sample] = np.squeeze(samples_mutations_load.loc[[sample], :].values)

for sample in female_LUSC_samples:
	samples_mutations_female_LUSC.loc[sample] = np.squeeze(samples_mutations_load.loc[[sample], :].values)



male_LUAD_cosine = DistanceMatrix(cosine_distances(samples_mutations_male_LUAD.values), samples_mutations_male_LUAD.index)
male_LUAD_tree_cosine = nj(male_LUAD_cosine)

print('Tree for males (LUAD) based on Cosine distance:')
print(male_LUAD_tree_cosine.ascii_art())
print('\n\n')

female_LUAD_cosine = DistanceMatrix(cosine_distances(samples_mutations_female_LUAD.values), samples_mutations_female_LUAD.index)
female_LUAD_tree_cosine = nj(female_LUAD_cosine)

print('Tree for females (LUAD) based on Cosine distance:')
print(female_LUAD_tree_cosine.ascii_art())
print('\n\n')

male_LUSC_cosine = DistanceMatrix(cosine_distances(samples_mutations_male_LUSC.values), samples_mutations_male_LUSC.index)
male_LUSC_tree_cosine = nj(male_LUSC_cosine)

print('Tree for males (LUSC) based on Cosine distance:')
print(male_LUSC_tree_cosine.ascii_art())
print('\n\n')

female_LUSC_cosine = DistanceMatrix(cosine_distances(samples_mutations_female_LUSC.values), samples_mutations_female_LUSC.index)
female_LUSC_tree_cosine = nj(female_LUSC_cosine)

print('Tree for females (LUSC) based on Cosine distance:')
print(female_LUSC_tree_cosine.ascii_art())
print('\n\n')



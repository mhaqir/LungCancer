

import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from glob import glob
import numpy as np
# import os

GRAPH_DIR = '/data/haghirebrahimm2/CTPsingle_AA_Lung/CONETT.txt'
CTP_DIR = '/data/haghirebrahimm2/CTPsingle_AA_Lung'

# clinical = pd.read_excel('Clinical.xlsx', 'Clinical Data')
# samples_excel = clinical['Tumor_Sample_Barcode'].values


Samples_after_CTP = [name[name[0:-1].rfind('/') + 1: -1] for name in glob(CTP_DIR + '/*/')]

with open(GRAPH_DIR, 'r') as f:
	lines = f.readlines()
lines = [line.rstrip().split(' ') for line in lines[1:]]
CONETT_graph = np.asarray(lines)


# Some samples absent from CONETT graph because CTPsingle assigns all of their mutations to one cluster
patients = list(set([item[0] for item in lines]))
diff_p = [p for p in Samples_after_CTP if p not in patients]



# genes1 = list(set([item[1] for item in lines]))
# genes2 = list(set([item[3] for item in lines]))
# diff_genes = [gene for gene in genes1 if gene not in genes2] + [gene for gene in genes2 if gene not in genes1]
# print(len(diff_genes))
# print(diff_genes)


# var_types1 = set(list(CONETT_graph[:, 2]))
# print(var_types1)
# var_types2 = set(list(CONETT_graph[:, 4]))
# print(var_types2)


# Checking for transitivity
for patient in patients:
	print('{}-----------------------------------------------------'.format(patient))
	sub_graph = CONETT_graph[CONETT_graph[:,0] == patient]
	graph_dict = {}
	for item in list(sub_graph):
		if (item[1], item[2], item[3]) not in graph_dict.keys():
			graph_dict[(item[1], item[2], item[3])] = [(item[3], item[4], item[5])]
		else:
			graph_dict[(item[1], item[2], item[3])] = graph_dict[(item[1], item[2], item[3])] + [(item[3], item[4], item[5])]

	first_second = list(zip(list(sub_graph[:,1]), list(sub_graph[:,2]), list(sub_graph[:,3])))
	for i in range(np.shape(sub_graph)[0]):
		if (sub_graph[i, 3], sub_graph[i, 4], sub_graph[i, 5]) in first_second:
			indices = [k for k, x in enumerate(first_second) if x == (sub_graph[i, 3], sub_graph[i, 4], sub_graph[i, 5])]

			for idx in indices:			
				if (sub_graph[idx, 3], sub_graph[idx, 4], sub_graph[idx, 5]) not in graph_dict[(sub_graph[i, 1], sub_graph[i, 2], sub_graph[i, 3])]:
					print(sub_graph[i,:])
					print(sub_graph[idx,:])
					print('------------------------')


# CONETT input graph
# edges = [(item[1], item[3]) for item in lines if item[0] == 'combined_SC284909']
# G = nx.DiGraph()
# G.add_edges_from(edges)
# pos = nx.random_layout(G)
# plt.figure(figsize=(15,15))
# nx.draw(G, pos, with_labels=True,node_size=1000,font_size=24) #
# plt.savefig('SC284909_graph.png')

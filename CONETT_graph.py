

import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from glob import glob
import numpy as np
# import os

GRAPH_DIR = '/data/haghirebrahimm2/CTPsingle_AA_Lung/CONETT_old.txt'
CTP_DIR = '/data/haghirebrahimm2/CTPsingle_AA_Lung'

# clinical = pd.read_excel('Clinical.xlsx', 'Clinical Data')
# samples_excel = clinical['Tumor_Sample_Barcode'].values


Samples_after_CTP = [name[name[0:-1].rfind('/') + 1: -1] for name in glob(CTP_DIR + '/*/')]

with open(GRAPH_DIR, 'r') as f:
	lines = f.readlines()
lines = [line.rstrip().split(' ') for line in lines[1:]]
CONETT_graph = np.asarray(lines)


# samples absent from CONETT graph because CTPsingle assigns all of their mutations to one cluster
patients = list(set([item[0] for item in lines]))
diff_p = [p for p in Samples_after_CTP if p not in patients]
# print(diff_p)
# print(len(diff_p))


# genes1 = list(set([item[1] for item in lines]))
# genes2 = list(set([item[3] for item in lines]))
# diff_genes = [gene for gene in genes1 if gene not in genes2] + [gene for gene in genes2 if gene not in genes1]
# print(len(diff_genes))
# print(diff_genes)

# CONETT inout graph
# edges = [(item[1], item[3]) for item in lines]
# G = nx.Graph()
# G.add_edges_from(edges)
# plt.figure(figsize=(25,25))
# nx.draw(G,node_size=60,font_size=8)
# plt.savefig('CONETT_graph.png')




# print(len(graph_dict))
# print(len(lines))
# print(list(set(graph_dict['ARID2'])))


for patient in patients:
	sub_graph = CONETT_graph[CONETT_graph[:,0] == patient]
	graph_dict = {}
	for item in list(sub_graph):
		if item[1] not in graph_dict.keys():
			graph_dict[item[1]] = [item[3]]
		else:
			graph_dict[item[1]] = graph_dict[item[1]] + [item[3]]
	for i in range(np.shape(sub_graph)[0]):
		if sub_graph[i, 3] in list(sub_graph[:,1]):
			if (sub_graph[i, 3] not in graph_dict[sub_graph[i, 1]]) and \
			(sub_graph[i, 1] not in graph_dict[sub_graph[i, 3]]):
				print(sub_graph[i,:])


				# if sub_graph[i, 0] not in non_transitive:
				# 	non_transitive.append(sub_graph[i, 0])



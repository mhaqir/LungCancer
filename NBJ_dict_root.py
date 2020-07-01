#!/bin/python3


from glob import glob
import sys
import numpy as np
import pandas as pd
from time import time
import itertools

from sklearn.metrics.pairwise import cosine_distances
from sklearn.neighbors import DistanceMetric


# from skbio import DistanceMatrix
from skbio.tree import nj
from skbio import DistanceMatrix

from Bio import Phylo
from ete3 import Tree

from Bio.Phylo.PhyloXML import Phylogeny
# from Bio.Phylo.TreeConstruction import DistanceMatrix
import networkx
from networkx.drawing import nx_agraph
networkx.graphviz_layout = nx_agraph.graphviz_layout


from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
# import Bio.Phylo.BaseTree.Tree as BTree
import pylab
from io import StringIO



# from hcluster import linkage, to_tree

clinical = pd.read_excel('Clinical.xlsx', 'Clinical Data')
# samples_excel = list(clinical['Tumor_Sample_Barcode'].values)
# samples_gender = list(clinical['Gender'].values)
# samples_histology = list(clinical['Histology'].values)

samples_id_gen_hist = clinical[['Tumor_Sample_Barcode', 'Gender', 'Histology']].values
# print(samples_id_gen_hist)

store = pd.HDFStore('samples_mutations.h5')
samples_mutations_load = store['samples_mutations']
samples_mutations_load_dict = samples_mutations_load.T.to_dict('list')

samples_mutations = {}
for k, v in samples_mutations_load_dict.items():
	new_k = k + '_' + samples_id_gen_hist[samples_id_gen_hist[:,0] == k, 1][0] + '_' + samples_id_gen_hist[samples_id_gen_hist[:,0] == k, 2][0]
	samples_mutations[new_k] = v



# print(sum(samples_mutations_load_dict['combined_SC284858']))
# print(sum(samples_mutations['combined_SC284858_Male_BAC']))


samples_mutations_values = np.asarray([v for k, v in samples_mutations.items()])
root_node = np.zeros((1, np.shape(samples_mutations_values)[1]), dtype = np.int8)
samples_mutations_values = np.concatenate((samples_mutations_values, root_node), axis = 0)
# print(np.shape(samples_mutations_values))
samples_mutations_keys = [k for k ,v in samples_mutations.items()]
samples_mutations_keys.append('root_node')
# print(samples_mutations_keys)
# print(type(samples_mutations_keys))
# print(type(samples_mutations_keys[0]))
# print(len(samples_mutations_keys))
# print(np.shape(samples_mutations_values))
# print(np.shape(samples_mutations_values))
# print(len(samples_mutations_keys))


all_cosine = DistanceMatrix(cosine_distances(samples_mutations_values), samples_mutations_keys)
# all_cosine_dist = cosine_distances(samples_mutations_values)

# dist_df = pd.DataFrame(all_cosine_dist, columns = samples_mutations_keys, index = samples_mutations_keys)

all_tree_cosine = nj(all_cosine)
all_tree_cosine.rooted = True
print(all_tree_cosine.root())
# print(all_tree_cosine.root_at('root_node'))

# print('Tree for all based on Cosine distance:')
# print(all_tree_cosine.ascii_art())

# print(all_cosine)

# all_tree_cosine = nj(all_cosine)


# print(type(all_tree_cosine))
# all_tree_cosine.write('all_tree.nwk', format ='newick')

# all_tree_cosine = all_tree_cosine[0:-1] + 'root;'
# print(all_tree_cosine)
# print(all_tree_cosine_newick)

# t = Tree(all_tree_cosine, format = 1)
# t.write(format = 1, outfile = 'all_tree.nwk')

# print(t)

# Phylo.write(t, 'all_tree.nwk', 'newick')
# Phylo.convert('all_tree.nwk', 'newick', 'all_tree.xml', 'nexml')


# tree = Phylo.read('all_tree.nwk', 'newick')

# constructor = DistanceTreeConstructor()
# upgmatree = constructor.upgma(all_cosine)

# print(upgmatree)

# Phylo.draw_ascii(tree)

# Phylo.draw.graphviz(tree, prog = 'dot')

# pos = list(zip(list(range(0, np.shape(all_cosine_dist)[0])), list(range(1, np.shape(all_cosine_dist)[1] + 1))))

# print(pos)
# dist = [list(all_cosine_dist[p[0], 0:p[1]]) for p in pos]


# print(type(samples_mutations_keys))


# from Bio.Phylo.TreeConstruction import _Matrix


# m = _Matrix(samples_mutations_keys, dist)

# dm = DistanceMatrix(samples_mutations_keys)
# print(dm)
# for s1, s2 in itertools.combinations(samples_mutations_keys, 2):
	# dm[s1, s2] = dist_df.loc[s1, s2]

# upgmatree = constructor.upgma(dm)

# print(upgmatree)
# Phylo.draw_ascii(upgmatree)
# Phylo.draw.graphviz(upgmatree, prog = 'dot')

# t = Tree(upgmatree, format = 1)
# t.write(format = 1, outfile = 'all_tree.nwk')

# njtree = constructor.nj(dm)
# njtree.rooted = True
# Phylo.write(njtree, 'njtree.xml', 'phyloxml')
# Phylo.convert('njtree.xml', 'phyloxml', 'njtree.nwk', 'newick')
# njtree_nwk = Phylo.read('njtree.nwk', 'newick')
# print(njtree_nwk)
# print(type(njtree))
# Phylo.draw_ascii(njtree)

# Phylo.draw(njtree, branch_labels = lambda c: c.branch_length)

# Phylo.draw_graphviz(njtree)
# pylab.show()
# pylab.savefig('tree.png')


# G = Phylo.to_networkx(njtree)
# print(G.node)

# g2 = networkx.convert_node_labels_to_integers(G, label_attribute='old_label' )

# Phylo.draw_graphviz(G)
# networkx.draw(G)
# pylab.show()
# pylab.savefig('tree1.png')


# t = Tree(njtree)

# phyloxml_tree = Phylo.BaseTree.Tree.as_phyloxml(njtree)
# print(type(phyloxml_tree))

# t = Tree(phyloxml_tree)

# Phylo.write(phyloxml_tree, 'tree.xml', 'phyloxml')

# Phylo.convert('tree.xml', 'phyloxml', 'njtree.nwk', 'newick')

# t = Phylo.read('njtree.nwk', 'newick')

# t = Tree(StringIO(t))

# tree = Phylogeny.from_tree(phyloxml_tree)

# for c in tree.get_terminals():
# 	# if c.name[0:5] == 'Inner':
# 	print(c.name)

# 	tree.c.color = 'blue'


# ml = tree.common_ancestor({'name': '*Male_LUAD'})
# ml.color = 'blue'



# print(tree)





import sys
from ete3 import Tree

t = Tree('{input.treefile}', format=0)

# Rename leaf names by removing '_1'
for leaf in t.iter_leaves():
    if leaf.name.endswith('_1'):
        leaf.name = leaf.name[:-2]

t.render('{output.image}', w=800, units='px')
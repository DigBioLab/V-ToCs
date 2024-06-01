# -*- coding: utf-8 -*-
"""
Created on Sun May 17 19:31:31 2020

@author: kosta
"""

import pandas as pd
import itertools
from Bio import Phylo
import os


def Newick_to_matrix():
    path = 'data'
    
    t = Phylo.read(f'{path}/Seqtreeupgma.txt', 'newick')
    
    d = {}
    for x, y in itertools.combinations(t.get_terminals(), 2):
        v = t.distance(x, y)
        d[x.name] = d.get(x.name, {})
        d[x.name][y.name] = v
        d[y.name] = d.get(y.name, {})
        d[y.name][x.name] = v
    for x in t.get_terminals():
        d[x.name][x.name] = 0
    
    m = pd.DataFrame(d)
    m.to_csv(f'{path}/Dataframe_Seq_upgma.csv', sep='\t')
    
if __name__ == "__main__":
    Newick_to_matrix()
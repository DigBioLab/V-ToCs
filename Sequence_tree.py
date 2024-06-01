# -*- coding: utf-8 -*-
"""
Created on Fri May 15 15:34:54 2020

@author: kosta
"""

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import os


def Sequence_tree():
    
    path = 'data'
    
    handle = open(f'{path}/aln-phylip.txt','r')
    aln = AlignIO.read(handle,'phylip')
    handle.close()
    
    calculator = DistanceCalculator('blosum62')
    dm = calculator.get_distance(aln)
    
    
    constructor = DistanceTreeConstructor()
    treenj = constructor.nj(dm)
    treeupgma = constructor.upgma(dm)
    
    Phylo.write(treenj, f'{path}/Seqtreenj.txt', 'newick')
    Phylo.write(treeupgma, f'{path}/Seqtreeupgma.txt', 'newick')



if __name__ == "__main__":
    Sequence_tree()
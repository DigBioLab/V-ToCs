# -*- coding: utf-8 -*-
"""
Created on Sun May 17 19:36:40 2020

@author: kosta
"""

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.spatial.distance as ssd
from scipy.cluster import hierarchy
from sklearn.cluster import SpectralClustering



#https://github.com/scipy/scipy/issues/8274
def getNewick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick



def Matrix_to_newick():
    #check log
    # =============================================================================
    # fold = r'C:\Users\kosta\Desktop\DTU\PhD\Projects\Venom_targetID\Comparison'
    # path = os.path.join(fold, 'pylog.txt')
    #
    # c = 0
    # with open(path, 'r') as fh:
    #     with open('log.txt','w') as fout:
    #         for line in fh:
    #             fout.write(fh.readline())
    #             c+=1
    #             if c== 500:
    #                 break
    # =============================================================================
    
    # save original working dir
    orig_dir = os.getcwd()
    os.chdir('data')
    
    file = 'Dataframe_TM.csv'
    file2 = 'Dataframe_RMS.csv'
    file3 = 'Dataframe_Seq_upgma.csv'
    
    df_tm = pd.read_csv(file, index_col = 0, sep='\t')
    df_rms = pd.read_csv(file2, index_col = 0, sep='\t')
    df_seq = pd.read_csv(file3, index_col = 0, sep='\t')
    
    if not os.path.exists('Figures'):
        os.mkdir('Figures')
    
    # =============================================================================
     #heatmap, TMalign score
    sns.heatmap(df_tm, xticklabels=False, yticklabels=False)
    plt.savefig('Figures/Heat_TMalign.png', dpi =300, format='png')
    #plt.savefig('Figures/Heat_TMalign.svg', format='svg')
    #
    sns.heatmap(df_rms, xticklabels=False, yticklabels=False, cmap='OrRd')
    plt.savefig('Figures/Heat_RMSD.png', dpi =300, format='png')
    #plt.savefig('Figures/Heat_RMSD.svg', format='svg')
    #
     #heatmap, RMSD score
    g = sns.clustermap(df_tm, col_cluster = False, xticklabels=False, yticklabels=True,
                       cbar_pos=(.2, .85, .2, .03), cbar_kws={"orientation": "horizontal"})
    g.ax_heatmap.tick_params(axis='both', which='both', length=0)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 0.5)
    g.fig.set_size_inches(42,40)
    g.savefig('Figures/Cluster_TMalign.png', dpi =200, format='png')
    link = g.dendrogram_row.linkage
    
    #
    # #clustermap, RMSD score
    g = sns.clustermap(df_rms, col_cluster = False, xticklabels=False, yticklabels=True, cmap='OrRd',
                       cbar_pos=(.2, .85, .2, .03), cbar_kws={"orientation": "horizontal"})
    g.ax_heatmap.tick_params(axis='both', which='both', length=0)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 0.5)
    g.fig.set_size_inches(42,40)
    g.savefig('Figures/Cluster_RMSD.png', dpi =200, format='png')
    link = g.dendrogram_row.linkage
    # =============================================================================
    
    
        
    vec_rms = ssd.squareform(df_rms)
    lrms = hierarchy.linkage(vec_rms, method='average')
    
    df_newtm = df_tm.copy()
    for i in df_newtm.columns:
        df_newtm[i][i] = 0
    print(df_newtm)
    vec_tm = ssd.squareform(df_newtm)
    ltm = hierarchy.linkage(1- vec_tm, method='average')
    
    tree_rms = hierarchy.to_tree(lrms)
    tree_tm = hierarchy.to_tree(ltm)
    
    rms = getNewick(tree_rms, '', tree_rms.dist, list(df_rms.columns))
    tm = getNewick(tree_tm, '', tree_tm.dist, list(df_newtm.columns))
    
    with open('RMSD_newick.txt', 'w') as fh:
        fh.write(rms)
    with open('TMalign_newick.txt', 'w') as fh:
        fh.write(tm)  
        
    # Return to the main working directory
    os.chdir(orig_dir)
    
    
if __name__ == "__main__":
    Matrix_to_newick()


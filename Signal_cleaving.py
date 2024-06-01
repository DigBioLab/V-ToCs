# -*- coding: utf-8 -*-
"""

"""

# pdb-tools
#!pip install pdb-tools

import os
import subprocess
import shutil
# import local script
from read_gff import read_gff


def cleave_pdb(protein_path, mod_path, res):
    # removing signal peptide residues
    subprocess.run(["pdb_delres", f"-:{res}", protein_path], 
                   stdout=open(mod_path, "w"), text=True)
    # renumbering residues
    subprocess.run(["pdb_reatom", mod_path], 
                   stdout=open(f"{mod_path}_temp", "w"), text=True)
    # renumbering atoms
    subprocess.run(["pdb_reres", f"{mod_path}_temp"], 
                   stdout=open(mod_path, "w"), text=True)
    # removing the temporary file
    os.remove(f"{mod_path}_temp")
    

def Signal_cleaving():
    path = 'data'
    pdb_path = os.path.join(path, 'Pdb_files')
    # filename for a simpler output list of proteins and respective signal cleavage sites
    cleave_list_file = os.path.join(path,'signal_cleaver_result.csv')

    # import cleavage site data from gff3
    cleave_list = read_gff(path, cleave_list_file)

    # directory for processed structures
    os.makedirs(f'{pdb_path}/modified', exist_ok=True)

    # parse the list
    for row_num, [protein_id, res] in cleave_list.iterrows():
        
        protein_path = os.path.join(pdb_path, f'{protein_id}_alphafold.pdb')
        mod_path = f'{pdb_path}/modified/{protein_id}.pdb'

        # cleave the residues
        if os.path.isfile(protein_path):
            cleave_pdb(protein_path, mod_path, res)
            
            # follow progress
            processed = len(os.listdir(f'{pdb_path}/modified'))
            total = len(os.listdir(f'{pdb_path}')) - 1
            if processed % 25 == 0:
                print(f'{processed/total*100:.0f}%')
              
       
    # copy the rest of proteins

    for file in os.listdir(pdb_path):
        protein_path = os.path.join(pdb_path, file)
        if os.path.isfile(protein_path):    # make sure the path doesn't refer to a folder
            protein_id, _ = file.split('_')
            mod_path = os.path.join(pdb_path, f'modified/{protein_id}.pdb')
            
            # only copy the structures that aren't in the folder
            if not os.path.isfile(mod_path):
                shutil.copy2(protein_path, mod_path)
    
    processed = len(os.listdir(f'{pdb_path}/modified'))
    total = len(os.listdir(f'{pdb_path}')) - 1
    print(f'{processed/total*100:.0f}%')
            
            
      
if __name__ == '__main__':
    Signal_cleaving()
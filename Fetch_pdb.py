# -*- coding: utf-8 -*-
"""


"""

import os
import requests


def Fetch_pdb():
    orig_dir = os.getcwd()
    os.chdir('data')
    fold = 'Uniprot_files'
    total_files = len(os.listdir(fold))
    
    count = 0
    listID = os.listdir(fold)
    no_struct = open('Unknown_structures.txt', 'w')         # Alphafold not available
    
    fold = r'Pdb_files'
    fold = os.path.normpath(fold)
    if not os.path.exists(fold):
        os.mkdir(fold)
    os.chdir(fold)
    
    for file in listID:
        url = r'https://alphafold.ebi.ac.uk/files/AF-' + file[:-4] + '-F1-model_v2.pdb'
        r = requests.get(url)
        name = file[:-4] + '_alphafold.pdb'
        if r.content.startswith(b'HEADER'):
            open(name, 'wb').write(r.content)
        else:
            no_struct.write(file[:-4] + '\n')
        count += 1
        print(f'{count / total_files * 100:.2f}%')
    no_struct.close()

    # Return to the main working directory
    os.chdir(orig_dir)

if __name__ == "__main__":
    Fetch_pdb()

#!/usr/bin/env python3

# -*- coding: utf-8 -*-


import os
import requests
import re


def Fetch_uniprot():
    # save original working dir
    orig_dir = os.getcwd()
    
    # create the 'data' directory
    if not os.path.exists('data'):
        os.mkdir('data')
    os.chdir('data')

    # download both Uniprot format and fasta files
    for format in ['txt', 'fasta']:
        url = f"https://rest.uniprot.org/uniprotkb/stream?format={format}" \
            "&query=%28%28taxonomy_id%3A8570%29%20AND%20%28reviewed%3Atrue%29%20" \
            "AND%20%28cc_tissue_specificity%3Avenom%29%29"
        r = requests.get(url)
        open(f'uniprot.{format}', 'wb').write(r.content)

    # make a directory for all Uniprot files (if not present)
    if not os.path.exists('Uniprot_files'):
        os.mkdir('Uniprot_files')

    # save uniprot files separately
    with open('uniprot.txt', 'r') as f:
        AC = None               # Uniprot ID
        for line in f:
            if line.startswith('ID'):
                ID = line

            # AC is None for the case of 2 or more AC lines in the Uniprot file
            elif line.startswith('AC') and AC is None:
                ACregex = re.match(r'AC\s+(\w+);', line)
                if ACregex:
                    AC = ACregex.group(1)
                    f_out = open('Uniprot_files/' + AC + '.txt', 'w')
                    f_out.write(ID)
                    f_out.write(line)

            elif line == '//\n':
                f_out.write(line)
                f_out.close()
                AC = None
            else:
                f_out.write(line)
    
    # Return to the main working directory
    os.chdir(orig_dir)
    
    
def optimize_fasta():
    # save original working dir
    orig_dir = os.getcwd()
    os.chdir('data')
    
    with open('uniprot.fasta') as f:
        long_id_count = 0
        long_id_file = open('Long_Uniprot_IDs.txt', 'w')
        optimized_fasta_file = open('uniprot_edit.fasta', 'w')
    
        for line in f:
    
            if line.startswith('>'):
                regex_id = re.match(r'>sp\|(\w+)(\|.*)', line)
                uniprot_id = regex_id.group(1)
    
                if len(uniprot_id) > 10:
                    long_id_count += 1
                    new_id = f'ID{long_id_count:04d}'
                    print(new_id, uniprot_id, sep='\t', file=long_id_file)
    
                    new_id_line = '>' + new_id + '\n'
                    optimized_fasta_file.write(new_id_line)
    
                else:
                    new_id_line = '>' + uniprot_id + '\n'
                    optimized_fasta_file.write(new_id_line)
    
            else:
                optimized_fasta_file.write(line)
    
        long_id_file.close()
        optimized_fasta_file.close()

    # Return to the main working directory
    os.chdir(orig_dir)
    
    
    
if __name__ == "__main__":
    Fetch_uniprot()
    optimize_fasta()



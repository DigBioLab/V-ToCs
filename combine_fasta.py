# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 18:43:10 2023

@author: Asus
"""

def export_fasta(protein_dict, filename):
    with open(filename, 'w') as f:
        for key in sorted(protein_dict.keys()): # sort keys alphabetically
            value = protein_dict[key]
            f.write('>' + key + '\n' + value + '\n')



def combine_fasta():
    
    path = 'data'
    sequences = {}
    cleaved_seq = []
    

    with open(f'{path}/processed_entries.fasta', "r") as fasta_file:
        current_sequence = ""
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                current_sequence = line[1:]
                cleaved_seq.append(current_sequence)
                sequences[current_sequence] = ""
            else:
                sequences[current_sequence] += line
                
    with open(f'{path}/uniprot_edit.fasta', "r") as fasta_file:      
        current_sequence = ""
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                current_sequence = line[1:]
                if current_sequence not in cleaved_seq:
                    sequences[current_sequence] = ""
            else:
                if current_sequence not in cleaved_seq:
                    sequences[current_sequence] += line

    export_fasta(sequences, f'{path}/combined.fasta')

    
if __name__ == '__main__':
    combine_fasta()
    
    
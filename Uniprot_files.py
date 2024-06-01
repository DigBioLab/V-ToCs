#!/usr/bin/env python3

# -*- coding: utf-8 -*-


import os
import requests
import re


def Uniprot_files():

    orig_dir = os.getcwd()
    os.chdir('data')
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
    
if __name__ == "__main__":
    Uniprot_files()
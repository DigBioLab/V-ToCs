# VToCs - Streamlit version

## Colab link

[Open in Colab](https://colab.research.google.com/drive/1ya3DjE2tyZsWWaV6TlWEgTrnwUDwpqn4?usp=sharing)

## Local installation

Clone the GitHub repository:  
```
git clone https://github.com/DigBioLab/V-ToCs.git
cd V-ToCs
```

Unpack uniprot files and unzip all archives in `data` directory:
```
python Uniprot_files.py
unzip 'data/*.zip' -d data/
```

Create a new virtual environment with python 3.9:

```
conda create -n vtocs-env python=3.9
conda activate vtocs-env
```

Install all dependencies: 

```
pip install -r requirements.txt
```

Run the streamlit app:
```
streamlit run Tree_gen_streamlit.py
```

When running on a machine with no screen (e.g. Google Colab), use the `offscreen` argument:
```
streamlit run Tree_gen_streamlit.py offscreen
```

It should start in your browser.  
Results are saved in `V-ToCs/Results`.
 

## Update V-ToCs

In order to keep the V-ToCs tool up to date with all the snake venom toxins currently available on Uniprot, follow the steps below.
To use the scripts, they should all be in the same directory. Change the working directory to the folder containing the scripts to make sure they work correctly.

### 1.	Fetch_uniprot.py

```
python Fetch_uniprot.py
```

Downloads Uniprot text files for each reviewed snake venom toxin, and a common FASTA file for all SVTs. 
The script will create directory `data/Uniprot_files`, where the Uniprot text files will be saved with names `[uniprotID].txt`. The FASTA file will be saved in the folder `data` named `uniprot.fasta`.
It will also edit the FASTA file so it can be used for sequence alignment in Clustal Omega. The edited file will be saved as `uniprot_edit.fasta` in the same directory `data`. In case there are Uniprot IDs which are longer than 10 symbols (maximum acceptable length is 10 for Clustal Omega), the longer IDs are substituted by shorter ones, and the list of substituted IDs is saved in a dictionary as `Long_Uniprot_IDs.txt`

### 2.	SignalP-6.0 

Many proteins on Uniprot contain signal peptides in their available structures and sequences. We remove them before all calculations. Use SignalP-6.0 to predict the cleavage sites of signal peptides on the proteins that contain signal peptides. Use any of the two options below:
 1)	DTU server

•	Go to https://services.healthtech.dtu.dk/services/SignalP-6.0/

•	Upload the file `uniprot_edit.fasta`

•	Choose organism: Eukarya

•	Choose output format: Short output (only text)

•	Choose model mode: Slow (fast can also work, but slow doesn’t take that long)

•	Submit

2)	Biolib server:

•	Go to https://dtu.biolib.com/SignalP-6

•	Same settings as with DTU server

•	Output format: Text

Save the output. From the output, copy the files `output.gff3` (Processed entries gff3) and `processed_entries.fasta` (Processed entries fasta) in the directory `data`


### 3.	combine_fasta.py

```
python combine_fasta.py
```

Combines the sequences, that had a signal peptide removed, from the file `processed_entries.fasta` with the sequences, that did not contain a signal peptide, from the file `uniprot_edit.fasta` into one fasta file: `combined.fasta`

### 4.	Clustal Omega

Use Clustal Omega to generate a sequence alignment file in the PHYLIP format. The tool is available at https://www.ebi.ac.uk/Tools/msa/clustalo/
Upload the combined FASTA file `combined.fasta` and set the output format to PHYLIP. Save the result as `aln-phylip.txt` in the directory `data`.

### 5.	Sequence_tree.py

```
python Sequence_tree.py
```

Generates a Newick tree file from the sequence alignment file made by Clustal Omega.
The script will search for the file `data/aln-phylip.txt` created by Clustal Omega. It will save two Newick files `Seqtreenj.txt` and `Seqtreeupgma.txt` in the same directory `data`.

### 6.	Newick_to_matrix.py

```
python Newick_to_matrix.py
```

Generates a distance matrix file from the Newick file made by `Sequence_tree.py`.
Run this script and it will search for the file `data/Seqtreeupgma.txt`. It will save the matrix file as `data/Dataframe_Seq_upgma.csv`.

### 7.	Fetch_pdb.py

```
python Fetch_pdb.py
```

Downloads AlphaFold2 structure files for each SVT. The script will search for the Uniprot text files in the directory `data/Uniprot_files`. It will save the structure files as `[UniprotID]_alphafold.pdb` for AlphaFold2 structures in the directory `data/Pdb_files`.

### 8.	pdb-tools

```
pip install pdb-tools
```

Install pdb-tools. This package is required for the next script to work. Use this command:

### 9.	Signal_cleaving.py

```
python Signal_cleaving.py
```

Modifies the .pdb files by cleaving the signal peptides from the toxin structures if a signal peptide is present. 
Reads the file `data/output.gff3` to extract the list of proteins have to be cleaved and at which site. Removes the residues belonging to the signal peptides and renumbers the atoms and residues. Saves the modified .pdb files as `data/Pdb_files/modified/[UniprotID].pdb`. The structures that did not contain the signal peptides, are copied into the `Pdb_files/modified` directory, also named as `[UniprotID].pdb`

### 10.	Runalign_tm_rms.py

Runs the tool TMalign that calculates similarity between structures independently of their sequences, and generates distance matrix and RMSD matrix for the structures analyzed.
In order to run this script, the TMalign.exe file has to be present in the same directory as the script. TMalign is a C++ program and it might need to be compiled from the cpp file on your computer before you can use this script. More info on: https://zhanggroup.org/TM-align/
The script will search for the pdb files in the directory `data/Pdb_files`. It will save the TM score matrix as `Dataframe_TM.csv`, TM distance matrix as `Dataframe_TMdist.csv` and RMSD distance matrix as `Dataframe_RMS.csv`.
There are two versions of the script available: 

•	`Runalign_tm_rms.py` with serial execution. It may take several days to process all .pdb files.
`python Runalign_tm_rms.py`

•	`Runalign_tm_rms_multiproc.py` that uses the Python multiprocessing module ProcessPoolExecutor. This script works a lot faster, depending on the available CPU it can take several hours to run. It uses all CPU power available, so it is recommended to only use it if you know what you are doing.
```
python Runalign_tm_rms_multiproc.py
```

### 11.	Matrix_to_newick.py

```
python Matrix_to_newick.py
```

Generates Newick tree files from the structure distance matrices. 
The script will search for the files `data/Dataframe_TMdist.csv` and `data/Dataframe_RMS.csv`. It will save the Newick files as `TMalign_newick.txt` and `RMSD_newick.txt` respectively, in the directory `data`.

### 12.	Tree_gen.py

The script generates a png, pdf, and svg file of the pruned tree, as well as a text file outlining the uniprot IDs of all of the toxins included in the final tree, as well as which snakes they come from.
With the minimal input, the script uses the uniprot text files from the directory `data/Uniprot_files`, as well as the Newick tree file, the distance dataframe file and the FASTA file. By default, Tree_gen.py will use the sequence tree file `data/Seqtreenj.txt`, the distance dataframe file for sequence `data/Dataframe_Seq_upgma.csv` and the original FASTA file `data/combined.fasta`. The output mentioned above is saved in the directory `Results`.


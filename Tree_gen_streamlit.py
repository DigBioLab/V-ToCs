import os
import sys
if len(sys.argv) > 1:
    if sys.argv[1] == "offscreen":
        os.environ['QT_QPA_PLATFORM']='offscreen'               #https://github.com/dkoslicki/TAMPA/issues/12
    else:
        raise ValueError("Invalid argument. Use 'offscreen'")
import re
import sys
import time
import pandas as pd
import numpy as np
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, AttrFace, RectFace
from Bio import SeqIO
import seaborn as sns
import streamlit as st
import multiprocessing
import json
import shutil
import zipfile

def render_tree(tree_rms, out_folder, ts):
    tree_rms.render(os.path.join(out_folder,'Test_prune.png'), dpi = 400, tree_style=ts)
    tree_rms.render(os.path.join(out_folder,'Test_prune.svg'), tree_style=ts)
    tree_rms.render(os.path.join(out_folder,'Test_prune.pdf'), tree_style=ts) 
    return None


@st.cache_data
def generate_tree(
        fam_select,
        org_file,
        org_type_infile,
        tree_option,
        fasta_file,
        acc_main,
        dist_cut,
        org_type,
        org_annot,
        col_consistent,
        fam_annot,
        frag_annot,
        tree_file,
        distance_file
):
    st.write("")
    #'''initialization'''
    accs = []
    org_origin = []
    frag_check = {}
    out = None
    files = []
    classification = []

    #'''files used'''
    #folder containing all data files
    folder = 'data'
    #tree file
    tree_path = os.path.join(folder, tree_file)
    #distance file for homology
    distance_path = os.path.join(folder, distance_file)
    #folder containing uniprot files for family annotation
    folder_uniprot = os.path.join(folder, r'Uniprot_files')
    #fasta file, checking for fragments
    fasta_path = os.path.join(folder, fasta_file)
    
    
    #'''variables used'''
    #intrepro domains to check for in the uniprot files
    inter = ['IPR006586', 'IPR024079', 'IPR009003', 'IPR003571', 'IPR001304', 'IPR002223', 'IPR001211', 'IPR001762', 'IPR003572']
    #toxin families that the domains correspond to
    map_toxins = ['P-III metalloprotease', '', 'Serine Protease', '', 'C-type lectin', 'Kunitz-type inhibitor', '', 'Disintegrin', 'Cytotoxin']
    #classes of toxins used
    class_set = {'Unclassified', 'Unclassified 3FTx', 'Long-chain α-neurotoxin', 'Short-chain α-neurotoxin', 'Cytotoxin',
                'Kunitz-type inhibitor', 'P-I metalloprotease', 'P-II metalloprotease',
                'P-III metalloprotease', 'Disintegrin', 'PLA2', 'Elapid PLA2', 'Viperid PLA2',
                'C-type lectin', 'Serine Protease'}
    
    #if only 1 family selected fam = the family. If all are selected then famf = all and fam =False:
    if fam_select == class_set:
        fam = False
        fam_file = 'all'
        fam_lst = False
    elif len(fam_select) == 1:
        fam = list(fam_select)[0]
        fam_file = 'all'
        fam_lst = False
    elif fam_select:
        fam_file = False
        fam = False
        fam_lst = fam_select
    else:
        fam = False
        fam_file = 'all'
        fam_lst = False


    
    if fam == False and fam_file == 'all' and fam_lst == False and org_file == 'all' and acc_main == False:
        warning="Warning: While VToCs can plot all toxins in all organisms - the tree will be too large to display. It should still save correctly. It may also take a while to run"
        html_warning=f"""<div style="text-align: center;"><p style='background-color:#FFFCEB;
                                            color:#916C04;
                                            font-size:16px;
                                            border-radius:3px;
                                            line-height:60px;
                                            padding-left:17px;
                                            opacity:1'>
                                            {warning}</style>
                                            <br></p></div>""" 
        st.markdown(html_warning, unsafe_allow_html=True)
        display = False
    else:
        display = True

    #'''read files'''
    #fasta file
    records = SeqIO.parse(fasta_path, "fasta")
    #distance dataframe
    df = pd.read_csv(distance_path, sep='\t', index_col = 0)
    #cleanup of dataframe
    df.index = [i.split('_')[0] for i in df.index]
    df.columns = [i.split('_')[0] for i in df.columns]

    if 'TM' in distance_file:
        df = df.apply(lambda x: 1-x)

    #tree file
    if 'Seq' in tree_file:
        with open(tree_path, 'r') as fh:
            raw = fh.read()
            tree_rms = Tree(raw, format=3) #5 for rms and tm, 3 for sequence
    else:
        with open(tree_path, 'r') as fh:
            raw = fh.read()
            tree_rms = Tree(raw, format=5) #5 for rms and tm, 3 for sequence
    

    
    st.header("")


    if org_file != 'all':
        warning=f"Warning: When a list of organisms is provided, the results are annotated by the same type of organism. Current annotation: by {org_type_infile}"
        html_warning=f"""<div style="text-align: center;"><p style='background-color:#FFFCEB;
                                            color:#916C04;
                                            font-size:16px;
                                            border-radius:3px;
                                            line-height:60px;
                                            padding-left:17px;
                                            opacity:1'>
                                            {warning}</style>
                                            <br></p></div>"""
        st.markdown(html_warning,unsafe_allow_html=True)
        org_type = org_type_infile

    if org_type == 'species' and col_consistent:
        warning="Warning: The consistent color scheme is not supported when the results are annotated by species"
        html_warning=f"""<div style="text-align: center;"><p style='background-color:#FFFCEB;
                                            color:#916C04;
                                            font-size:16px;
                                            border-radius:3px;
                                            line-height:60px;
                                            padding-left:17px;
                                            opacity:1'>
                                            {warning}</style>
                                            <br></p></div>"""
        st.markdown(html_warning,unsafe_allow_html=True)
        col_consistent = False

    max_dist = df.max().max()

    os.makedirs('Results', exist_ok=True)
    timestamp = time.strftime('%Y_%m_%d_%H_%M_%S')
    out_folder = os.path.join('Results', timestamp)
    os.mkdir(out_folder)
    out_path = os.path.join(out_folder, 'Raw.txt')
    inputtxt_path = os.path.join(out_folder, 'input.txt')
    with open(inputtxt_path, 'w', encoding='utf-8') as fin:
        input_str = "".join([
            f"toxin families: fam_select = {fam_select}\n",
            f"organisms: org_file = {org_file}\n",
            f"type of organism: org_type_infile = {org_type_infile}\n",
            f"tree option: tree_option = {tree_option}\n",
            f"fasta file: fasta_file = {fasta_file}\n",
            f"uniprot id for homology distance: acc_main = {acc_main}\n",
            f"homology distance cutoff %: dist_cut = {dist_cut}\n",
            f"annotate by: org_type = {org_type}\n",
            f"Annotate organism by tree leaf colouring scheme?: org_annot = {org_annot}\n",
            f"Keep the genera colour scheme consistent?: col_consistent = {col_consistent}\n",
            f"Annotate toxin family by tree leaf colouring scheme?: fam_annot = {fam_annot}\n",
            f"Annotate fragments by a line next to the node?: frag_annot = {frag_annot}\n"
        ])
        print(input_str, file=fin)

        
    col1, col2, col3 = st.columns([1, 0.5, 1])

    with col2:
        with st.spinner('Wait for it...'):
            #tree with all orgs
            if org_file == 'all':
                for txtfile in os.listdir(folder_uniprot):
                    txtfile_path = os.path.join(folder_uniprot, txtfile)
                    with open(txtfile_path, 'r') as fin:
                        for line in fin:
                            if line.startswith('OS'):
                                species_line = line[5:-1]
                                break
            
                        if org_type == 'genus':
                            org = species_line.split()[0]
                        elif org_type == 'species':
                            org = ' '.join(species_line.split(' ')[:2])        # join and split to remove is subspecies name
                        else:
                            print("Invalid input for organism type. Options are either 'genus' or 'species'. Please try again")
                            sys.exit(2)
            
                        acc = txtfile[:-4]
                        accs.append(acc)
                        org_origin.append(org)
                orgs = list(set(org_origin))
            
            #tree with list of selected orgs
            else:
                orgs = org_file.splitlines()
                orgs = list(set(orgs))
                with open(out_path, 'w', encoding='utf-8') as fout:
            
                    for txtfile in os.listdir(folder_uniprot):
                        txtfile_path = os.path.join(folder_uniprot, txtfile)
                        with open(txtfile_path, 'r') as fin:
                            for line in fin:
                                if line.startswith('OS'):
                                    species_line = line[5:-1]
                                    break
            
                            if org_type == 'genus':
                                org = species_line.split()[0]
                            else:
                                org = ' '.join(species_line.split(' ')[:2])  # join and split to remove is subspecies name

            
                        if org in orgs:
                            acc = txtfile[:-4]
                            accs.append(acc)
                            org_origin.append(org)
                            print(org+'\t'+acc, file=fout)
            
            #cleanup of pdb name suffix in the tree node names
            for leaf in tree_rms.get_leaves():
                leaf.name = leaf.name.split('_')[0]
            
            #check for accessions in the tree, remove those that are not present
            for i in reversed(range(len(accs))):
                match = tree_rms.search_nodes(name=accs[i])
                if len(match) == 0:
                    accs.pop(i)
                    org_origin.pop(i)
            
            #prune the tree, keep only accessions of interest
            tree_rms.prune(accs)
            
            #check for fragments, annotate those that are fragements
            for rec in records:
                acc = rec.id
            
                if acc in accs:
                    st1 = rec.description
            
                    if 'Fragment' in st1:
                        frag_check[acc] = True
                    else:
                        frag_check[acc] = False

            
            #double check all accessions and whether they are in the fragment dictionary
            for acc in accs:
                if acc not in frag_check:
                    frag_check[acc] = False
            #initialize uniprot files to search for toxin family annotation
            files_to_search = [i+'.txt' for i in accs]
            
            #search files for toxin family annotation            
            
            for file in files_to_search:
                fh = open(os.path.join(folder_uniprot, file))
                content = fh.read()
                fh.close()

                file = file[:-4]

                for i, dom in enumerate(inter):
                    match = re.search(f'DR\s+InterPro;\s{dom};\s([\w\/-]+)\.', content)
                    
                    # identify through domain numbers
                    if match:
                        if dom == 'IPR024079':
                            match = re.search(f'DR\s+InterPro;\sIPR001762;\s([\w\/-]+)\.', content)
                            if not match:
                                out = 'P-I metalloprotease'

                            else:
                                match = re.search(f'DR\s+InterPro;\sIPR034027;\s([\w\/-]+)\.', content)
                                if match:
                                    out = 'P-III metalloprotease'
                                else:
                                    out = 'P-II metalloprotease'

                            classification.append(out)
                            break

                        elif dom == 'IPR001211':
                            match = re.search(r'OC\s{3}\S.+\.', content)
                            lst = match.group().split(';')
                            try:
                                if lst[2] == ' Viperidae':
                                    out = 'Viperid PLA2'
                                elif lst[2] == ' Elapidae':
                                    out = 'Elapid PLA2'
                                else:
                                    out = 'PLA2'
                            except:
                                print(file)
                            classification.append(out)
                            break

                        elif dom == 'IPR003571':
                            match = re.search(r'DR\s+InterPro;\sIPR003572;\s([\w\/-]+)\.', content)
                            if match:
                                out = 'Cytotoxin'

                            else:
                                match = re.search(r'DE\s{3}(RecName|AltName): Full=Long', content)
                                if match:
                                    out = 'Long-chain α-neurotoxin'
                                else:
                                    match = re.search(r'DE\s{3}(RecName|AltName): Full=Short', content)
                                    if match:
                                        out = 'Short-chain α-neurotoxin'
                                    else:
                                        out = 'Unclassified 3FTx'
                                        with open("unclas_3ftxs.txt", "a") as unclas_file:
                                            print(file, file=unclas_file)

                            classification.append(out)
                            break

                        else:
                            out = map_toxins[i]
                            classification.append(out)
                            break

                    # if dom didn't work, id through DE line
                    else:
                        if i == len(inter) - 1:
                            match_svmp = re.search(r'DE\s+Short=SVMP', content)
                            match_pla = re.search(r'DE\s+Short=(sv)?PLA2', content)
                            match_3ftx = re.search(r'DE\s{3}(RecName|AltName): Full=Three-finger', content)
                            match_ScNTx = re.search(r'DE\s{3}(RecName|AltName): Full=Short.+neurotoxin', content)
                            match_LcNTx = re.search(r'DE\s{3}(RecName|AltName): Full=Long.+neurotoxin', content)
                            match_cyto = re.search(r'DE\s{3}(RecName|AltName): Full=.*(C|c)ytotoxin', content)
                            if match_pla:
                                match = re.search(r'OC\s{3}\S.+\.', content)
                                lst = match.group().split(';')
                                if lst[2] == ' Viperidae':
                                    out = 'Viperid PLA2'
                                elif lst[2] == ' Elapidae':
                                    out = 'Elapid PLA2'
                                else:
                                    out = 'PLA2'
                            elif match_svmp:
                                out = 'P-II metalloprotease'        # is it P-II though?
                            elif match_ScNTx:
                                out = 'Short-chain α-neurotoxin'
                            elif match_LcNTx:
                                out = 'Long-chain α-neurotoxin'
                            elif match_cyto:
                                out = 'Cytotoxin'
                            elif match_3ftx:
                                out = 'Unclassified 3FTx'

                            else:
                                out = 'Unclassified'
                            classification.append(out)

                if out is None:
                    print(f'Not found, file {file}')
            
            # tree with selected toxin families
            
            if fam_file != 'all':
            
                indices = [i for i in range(len(classification)) if classification[i] in fam_lst]
                accs = [accs[i] for i in indices]
                org_origin = [org_origin[i] for i in indices]
                classification = [classification[i] for i in indices]
            
            # tree with one selected toxin family
            if fam:
            
                if fam not in classification:
                    print('Invalid family to filter for: Options are:')
                    for i in list(set(classification)):
                        print(i)
                    print('Please try again.')
                    sys.exit(5)
            
                indices = [i for i in range(len(classification)) if classification[i] == fam]
                accs = [accs[i] for i in indices]
                org_origin = [org_origin[i] for i in indices]
                classification = [classification[i] for i in indices]
            
            #prune the tree, keep only accessions of interest
            tree_rms.prune(accs)
            
            #distance homology calculations and styling
            if acc_main:
            
                if acc_main not in accs:
                    st.error('Accession to base homology on is not part of accessions included in the tree. Please try again with a valid accession or without homology calculations.')
                    st.stop()
            
                dist_lst = []
                #get homology distance values
                for acc in accs:
                    if acc != acc_main:
                        dist = df.loc[acc_main][acc]
                        if np.isnan(dist):
                            dist = df.loc[acc][acc_main]
                        hom = round(100 - ((dist/max_dist) *100),2)
                        if dist_cut:
                            dist_lst.append(hom)
                    else:
                        dist_lst.append(100)
            
                if dist_cut:
                        
                    indices = [i for i in range(len(dist_lst)) if dist_lst[i] > dist_cut]
                    accs = [accs[i] for i in indices]
                    org_origin = [org_origin[i] for i in indices]
                    classification = [classification[i] for i in indices]
                    dist_lst = [dist_lst[i] for i in indices]
            
            #prune tree based on distance
            tree_rms.prune(accs)

            #'''Tree style'''
            ts= TreeStyle()
            
            #colors used for tree annotation of organisms
            un_orgs = list(set(org_origin))
            if col_consistent == False:
                cmap = sns.color_palette('husl', len(un_orgs)) # seaborn create col palette with length of number of orgs
                cmap = cmap.as_hex()
            else:
                with open("data/genera_cmap.json") as f:
                    cmap = json.load(f)

            #colors used for tree annotation of families
            cmap_fam = {'Unclassified 3FTx':'#e3665b', 'Long-chain α-neurotoxin':'#995e59', 'Short-chain α-neurotoxin':'#e8c2be',
                        'Cytotoxin':'#941e12', 'PLA2':'#9c6fad', 'Viperid PLA2':'#58256b', 'Elapid PLA2':'#ddadf0',
                        'P-I metalloprotease':'#c6f0ad', 'P-II metalloprotease':'#76ad55', 'P-III metalloprotease':'#326315',
                        'C-type lectin':'#5f6ed9', 'Kunitz-type inhibitor':'#f7a257', 'Disintegrin':'#279c79',
                        'Serine Protease':'#86e3e3', 'Unclassified':'#cbd0d1'}
            
            #organism legend
            if org_annot:
                for k, org in enumerate(un_orgs):
                    face_leg = TextFace(org, ftype='Arial', fsize=20)
                    face_leg.margin_top = 10
                    face_leg.margin_right = 10
                    face_leg.margin_left = 10
                    face_leg.margin_bottom = 10
                    face_leg.opacity = 0.5 # from 0 to 1]
                    if col_consistent == False:
                        face_leg.background.color = cmap[k]
                    else:
                        face_leg.background.color = cmap[org]
                    ts.legend.add_face(face_leg, column=0)
            
            #toxin family legend
            if fam_annot:
                for fam in list(set(classification)):
            
                    face = TextFace(fam, ftype='Arial', fsize=20)
                    face.margin_top = 5
                    face.margin_right = 5
                    face.margin_left = 5
                    face.margin_bottom = 5
                    face.opacity = 0.5
                    face.background.color = cmap_fam[fam]
                    ts.legend.add_face(face, column=1)
            
            #node styling
            for i, acc in enumerate(accs):
            
                org = org_origin[i]
                family = classification[i]
                org_index = un_orgs.index(org)
            
                ns = NodeStyle()
                if org_annot:
                    if col_consistent == False:
                        ns["fgcolor"] = cmap[org_index]
                    else:
                        ns["fgcolor"] = cmap[org]
                    ns["size"] = 15
                if fam_annot:
                    ns["bgcolor"] = cmap_fam[family]
            
                node = tree_rms & acc
            
                if frag_annot:
                    #check if fragment, annotate with a big ass black line
                    if frag_check[acc]:
            
                        face_frag = RectFace(1000, 10, 'black', 'black', label=None)
                        node.add_face(face_frag, column= 3, position = 'branch-right')
            
                node.set_style(ns)
            
            if acc_main:
                #distance homology calculations and styling
                for i, acc in enumerate(accs):
                        node = tree_rms & acc
                        node.add_features(distance = dist_lst[i])
            
                #annotate homology distance
                for acc in accs:
            
                    node = tree_rms & acc
            
                    if node.distance < 30:
                        face_dist = AttrFace('distance', fsize= 6, fgcolor="#1A6D09")
                    elif node.distance < 50:
                        face_dist = AttrFace('distance', fsize= 8, fgcolor="#56F800")
                    elif node.distance < 75:
                        face_dist = AttrFace('distance', fsize= 10, fgcolor="#F89200")
                    else:
                        face_dist = AttrFace('distance', fsize= 12, fgcolor="#F80B00")
            
                    node.add_face(face_dist, column = 0, position = 'branch-top')
            
            #remove inner node styling
            for node in tree_rms.traverse():
                ns_inner = NodeStyle()
                ns_inner["size"] = 0
            
                if not node.is_leaf():
                    node.set_style(ns_inner)
            
            #final assignments
            for node in tree_rms.traverse():
                node.dist = 1
            ts.legend_position = 1
            ts.show_leaf_name = True
            ts.scale = 80
            ts.branch_vertical_margin = 20
            
            #create tree figures
            process = multiprocessing.Process(target=render_tree, args=(tree_rms, out_folder, ts))
            process.start()
            process.join()

            #render_tree(tree_rms, out_folder, ts)
            #raw file with accs and orgs
            with open(out_path, 'w', encoding='utf-8') as fout:
                for acc, org, cl in zip(accs, org_origin, classification):
                    print(f'{acc}\t{org}\t{cl}', file=fout)

            zfile_path = os.path.join(out_folder, f"{timestamp}.zip")
            # add files into zipfile
            with zipfile.ZipFile(zfile_path, 'w') as zfile:
                for f in os.listdir(out_folder):
                    if ".zip" not in f:
                        f_path = os.path.join(out_folder, f)
                        zfile.write(f_path, arcname=f)

    return out_folder, timestamp, display



def main():
    st.set_page_config(
         page_title="Digital Biotechnology V-ToCs: Venom TOxin CluStering",
         page_icon=":snake:",
         layout="wide",
         initial_sidebar_state="expanded",
         menu_items={
             'Get Help': 'https://digital-biotechnology.com/',
             'About': "https://digital-biotechnology.com/"
         }
    )    
    st.title(':snake: Welcome to V-ToCs: Venom TOxin CluStering')
    st.write("Tree figure generation, pruning and homology calculation.")
    form = st.form(key="submit_form")   #, clear_on_submit = True)
    
    with form:
        st.subheader(':point_right:  Input options:')
        fam_select = st.multiselect(
            #'Select which toxin families to prune with.',
            'Select which toxin families to prune with.',
            [
                "Unclassified", 
                "Unclassified 3FTx", 
                "Long-chain α-neurotoxin", 
                "Short-chain α-neurotoxin", 
                "Cytotoxin", 
                "Kunitz-type inhibitor", 
                "P-I metalloprotease", 
                "P-II metalloprotease", 
                "P-III metalloprotease", 
                "Disintegrin", 
                "PLA2", 
                "Elapid PLA2", 
                "Viperid PLA2", 
                "C-type lectin", 
                "Serine Protease"
            ]
        )
        fam_select = set(fam_select)

    
        st.markdown("---") 
        st.subheader(':evergreen_tree:  Tree options:')

        col1, col2 = st.columns(2, gap="large")
        with col1:

            org_type_infile = st.radio(
                'Plotting a list of organisms by:',
                ['genus', 'species'],
            )

            with open("data/genera_cmap.json", "r") as f:
                genera_list = list(json.load(f).keys())
            org_genera_file = st.multiselect(
                'Select which genera you would like to plot. Leave blank to enter them manually below.',
                genera_list
            )
            org_genera_file = set(org_genera_file)

            org_file_textarea_input = st.text_area('Enter the list of species or genera you would like to plot. Paste in one name per line. Leave blank to plot all available species/genera.')
            # remove unnecessary characters, capitalize, don't count empty lines
            org_file = []
            for el in org_file_textarea_input.split("\n"):
                org = el.strip().strip(",;")
                if org:
                    org = [s for s in org.split(" ") if s]
                    org = " ".join(org)
                    org = org.capitalize()
                    org_file.append(org)

            
            if org_type_infile == 'genus' and org_genera_file:
                org_file = "\n".join(org_genera_file)
            elif org_file:
                org_file = "\n".join(org_file)
            else:
                org_file = 'all'

            
        with col2:
            #fasta_file = st.file_uploader('Fasta file for fragment check. Upload your own file or leave blank to use the default all toxins fasta file.')
            #if not fasta_file:
            #    fasta_file = 'combined.fasta'
            fasta_file = 'combined.fasta'
            
            tree_option = st.selectbox(
                  'Select the tree file to work with:',
                  ('Sequence_blosum45','Sequence_blosum62','Sequence_pam250', 'Structure-TM', 'Structure-RMSD')
            )
            tree_and_distance_files  = {
                'Sequence_blosum62': ('Seqtreenj.txt','Dataframe_Seq_upgma.csv'), 
                'Sequence_blosum45': ('Seqtreenj_blosum45.txt','Dataframe_Seq_upgma_blosum45.csv'), 
                'Sequence_pam250': ('Seqtreenj_pam250.txt','Dataframe_Seq_upgma_pam250.csv'), 
                'Structure-TM': ('TMalign_newick.txt', 'Dataframe_TMdist.csv'), 
                'Structure-RMSD': ('RMSD_newick.txt','Dataframe_RMS.csv')
            }
            tree_file, distance_file = tree_and_distance_files.get(tree_option)

            acc_main = st.text_input('Get homology distance values based on input accession. Must be a uniprot accession  number.')
            acc_main = False if acc_main == '' else acc_main
            
            dist_cut = st.slider('Distance cutoff to prune tree with. Needs to be used in conjuction with homology distance based on input accession. This is a percentage number from 0-100.', 0, 100, 0)
            dist_cut = None if dist_cut == 0 else float(dist_cut)
            
        #Default organism type is genus, so you can use a genera list without specifying the organism type.
    
        st.markdown("---") 
        st.subheader(':pushpin:  Annotation options:')
        col1, col2 = st.columns(2, gap="large")
        with col1:
    
            org_type = st.radio(
                'Annotate the results by genus or species?',
                ['genus', 'species'],
                key=2,
                )
            
            org_annot = st.radio(
                 'Annotate organism genus or species by tree leaf colouring scheme?',
                 ('Annotate organism', 'No annotation'),
                 key=3,
                 )
            org_annot = True if org_annot == 'Annotate organism' else False
            
            col_consistent = st.radio(
                 #'Optimise the best colours for each individual run or keep them consistent? (Second option is useful if you want to run the tool many times and have consistent colouring across different graphs.)',
                 'Keep the genera colour scheme consistent?',
                 ('Keep consistent', 'No'),
                 key=4,
                 )
            col_consistent = True if col_consistent == 'Keep consistent' else False
    
        with col2:
            fam_annot = st.radio(
                 'Annotate toxin family by tree leaf colouring scheme?',
                 ('Annotate toxin family', 'No annotation'),
                 key=5
                 )
            fam_annot = True if fam_annot == 'Annotate toxin family' else False
            
            frag_annot = st.radio(
                'Annotate fragments by a line next to the node?',
                ('Annotate with line', 'No annotation'),
                key=6
                )
            frag_annot = True if frag_annot == 'Annotate with line' else False
    
        submitted = st.form_submit_button(label="Submit")
    

    if submitted:
        out_folder, timestamp, display = generate_tree(
                                                fam_select,
                                                org_file,
                                                org_type_infile,
                                                tree_option,
                                                fasta_file,
                                                acc_main,
                                                dist_cut,
                                                org_type,
                                                org_annot,
                                                col_consistent,
                                                fam_annot,
                                                frag_annot,
                                                tree_file,
                                                distance_file
        )
        
        #col1, col2, col3 = st.columns([1, 0.3, 1])

        #with col2:
            #st.cache_data.clear()
        zfile_path = os.path.join(out_folder, f"{timestamp}.zip")
        with open(zfile_path, "rb") as file:
            st.download_button(
                label="Download",
                data=file,
                file_name=f"v-tocs_tree_{timestamp}.zip",
                #mime="image/png",
                #on_click=submit_form
            )

        txt=f"Damn it worked. Surprising... Results are located in {out_folder}"   
        st.success(txt, icon="✅")

        if display:
            image = out_folder + '/Test_prune.png'
            print(image)
            st.image(image, use_column_width='always')


if __name__ == "__main__":
    main()                

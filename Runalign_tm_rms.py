import subprocess
import os
import re
import pandas as pd
import numpy as np


def Runalign_tm_rms():
    path = 'data'
    # files to analyze
    folder = f'{path}/Pdb_files/modified'
    files = []
    for file in os.listdir(folder):
        if file.endswith('pdb'):
            files.append(file)
    
    files = sorted(files)
            
    df = pd.DataFrame(np.nan, columns=files, index=files)           # original TM scores
    df_tm_dist = pd.DataFrame(np.nan, columns=files, index=files)   # distance scores from TM
    df_rms = pd.DataFrame(np.nan, columns=files, index=files)       # RMSD scores
    
    print(len(df), len(df.columns))
    print(len(df_rms), len(df_rms.columns))
    
    c = 0
    prog = 0
    targets = files.copy()
    
    print('Initialized')
    
    for t in files:
        progt = 0
        target = targets[0]
        df.loc[target][target] = 1
        df_tm_dist.loc[target][target] = 0
        df_rms[target][target] = 0
    
        if len(targets) > 1:
                
            for i in targets[1:]:
                #print(f'Comparing {i}, {target}')
                command = f'./TMalign.exe {folder}/{i} {folder}/{target}'.split(' ')
                res = subprocess.run(command, stdout=subprocess.PIPE)
                
                if res.returncode == 0:
                    #print(f'Successful {i}, {target}')
                    if c == 0:
                        with open('pylog.txt','w') as fh:
                            out = res.stdout.decode()
                            fh.write(out)
                            fh.write('Iteration 1\n')
                            scores = re.findall(r'TM-score= (\d*\.?\d*)',out)
                            df.loc[i][target] = float(scores[1])
                            df.loc[target][i] = float(scores[1])
                            df_tm_dist.loc[i][target] = 1 - float(scores[1])
                            df_tm_dist.loc[target][i] = 1 - float(scores[1])
                            rmsd = re.findall(r'RMSD=\s+(\d*\.?\d*)',out)
                            df_rms.loc[i][target] = float(rmsd[0])
                            df_rms.loc[target][i] = float(rmsd[0])
    
                        c += 1
                        progt += 1
                    
                    else:
                        with open('pylog.txt','a') as fh:
                            out = res.stdout.decode()
                            fh.write(out)
                            fh.write('\n')
                            scores = re.findall(r'TM-score= (\d*\.?\d*)',out)
                            try:
                                df.loc[i][target] = float(scores[1])
                                df.loc[target][i] = float(scores[1])
                                df_tm_dist.loc[i][target] = 1 - float(scores[1])
                                df_tm_dist.loc[target][i] = 1 - float(scores[1])
                            except IndexError:
                                df.loc[i][target] = 'Failed_re'
                                df.loc[target][i] = 'Failed_re'
                                df_tm_dist.loc[i][target] = 'Failed_re'
                                df_tm_dist.loc[target][i] = 'Failed_re'
                                print(f'Failed {target} and {i}, TMalign')
                            
                            rmsd = re.findall(r'RMSD=\s+(\d*\.?\d*)',out)
                            try:
                                df_rms.loc[i][target] = float(rmsd[0])
                                df_rms.loc[target][i] = float(rmsd[0])
                            except IndexError:
                                df_rms.loc[i][target] = 'Failed_re'
                                df_rms.loc[target][i] = 'Failed_re'
                                print(f'Failed {target} and {i}, rmsd')
                        
                        progt += 1
                
                else:
                    print(f'Problem with {i}, {target}')
                    df.loc[i][target] = 'Failed2_rc'
                    df.loc[target][i] = 'Failed2_rc'
                    df_tm_dist.loc[i][target] = 'Failed2_rc'
                    df_tm_dist.loc[target][i] = 'Failed2_rc'
                    df_rms.loc[i][target] = 'Failed_rc'
                    df_rms.loc[target][i] = 'Failed_rc'
            
            tpop = targets.pop(0)
            #print(targets,tpop)
    
        else:
            print('Length one, finished')
    
        print(f'Compared {progt} with {target}')
        prog += 1
        print(f'{prog} out of {len(files)}')
    
    print('Cleaning names')
    
    df.columns = [i[:-17] if 'sel' in i else i[:-4] for i in df.columns]
    df.index = [i[:-17] if 'sel' in i else i[:-4] for i in df.index]
    
    df_tm_dist.columns = [i[:-17] if 'sel' in i else i[:-4] for i in df_tm_dist.columns]
    df_tm_dist.index = [i[:-17] if 'sel' in i else i[:-4] for i in df_tm_dist.index]
    
    df_rms.columns = [i[:-17] if 'sel' in i else i[:-4] for i in df_rms.columns]
    df_rms.index = [i[:-17] if 'sel' in i else i[:-4] for i in df_rms.index]
    
    print('Writing dataframes')
    df.to_csv(f'{path}/Dataframe_TM.csv', sep='\t')
    df_tm_dist.to_csv(f'{path}/Dataframe_TMdist.csv', sep='\t')
    df_rms.to_csv(f'{path}/Dataframe_RMS.csv', sep='\t')
    
    print('Finished', c)
    
if __name__ == "__main__":
    Runalign_tm_rms()

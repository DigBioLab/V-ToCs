import subprocess
import os
import re
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed



def execute_command(command_full, count, total_count):
    command, i_protein, target = command_full
    res = subprocess.run(command, stdout=subprocess.PIPE)
    if res.returncode == 0:
        try:
            out = res.stdout.decode()
            TMscore = re.findall(r'TM-score= (\d*\.?\d*)',out)[1]
            TMscore = float(TMscore)
            TMdist = 1 - TMscore
            rmsd = re.findall(r'RMSD=\s+(\d*\.?\d*)',out)[0]
        except IndexError:
            TMscore = 'Failed_re'
            TMdist =  'Failed_re'
            rmsd =    'Failed_re'
            print(f'Failed {target} and {i_protein}, TMalign')
            
        # follow progress    
        if count % 300 == 0:
            print(f'{count/total_count*100:.2f}%; TMscore = {TMscore}')
            with open('runalign_log_processpool.txt', 'a') as f:
                f.write(f'{count/total_count*100:.2f}%; TMscore = {TMscore}\n')
        return (TMscore, TMdist, rmsd, i_protein, target, out)
    else:
        raise ValueError(f'execute_command failed: i = {i_protein}, target = {target}')
        #return None



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
    
    for i_df in [df, df_tm_dist, df_rms]:
        print(len(i_df), len(i_df.columns))
    
    
    c = 0   # count
    prog = 0
    targets = files.copy()
    
    print('Initialized')
    
    command_list = []
    
    for t in files:
        target = targets[0]
        df.loc[target][target] = 1
        df_tm_dist.loc[target][target] = 0
        df_rms[target][target] = 0
    
        if len(targets) > 1:
                
            for i in targets[1:]:
                #print(f'Comparing {i}, {target}')
                command = f'./TMalign {folder}/{i} {folder}/{target}'.split(' ')
                # argument with the command and with the proteins involved
                command_full = [command, i, target]
                command_list.append(command_full)            
            tpop = targets.pop(0)
            #print(targets,tpop)
            
        else:
            print('Length one, finished')
        prog += 1
        print(f'{prog} out of {len(files)}')
        
    total_count = len(command_list)
    print(f'{total_count} comparisons')
    
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(execute_command, 
                                   command_full, count, total_count) for 
                                   (count, command_full) in enumerate(command_list)]

    
        for future in as_completed(futures):
            output = future.result()
            
            # save scores
            if output is not None:
                TMscore, TMdist, rmsd, i, target, out = output
                # save TMalign outputs
                '''with open('pylog.txt','a') as fh:
                    #fh.write(out)
                    #fh.write('\n')'''
                
                # Original TM scores
                df.loc[i][target] = float(TMscore)
                df.loc[target][i] = float(TMscore)
                
                # Distance scores from TM
                df_tm_dist.loc[i][target] = float(TMdist)
                df_tm_dist.loc[target][i] = float(TMdist)
                
                # RMSD scores
                df_rms.loc[i][target] = float(rmsd)
                df_rms.loc[target][i] = float(rmsd)

            else:
                raise ValueError(f'execute_command failed: i = {i}, target = {target}') 


    
    
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


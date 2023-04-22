import os
import math
import sys
import traceback
import re
import shutil
from datetime import datetime
import pandas as pd
import numpy as np
import argparse
import subprocess
import linecache

if __name__ == "__main__":
    seq_file_dire = sys.argv[1]
    subprocess.run('cat '+seq_file_dire+'*Information.txt > '+seq_file_dire+'all_Infor.txt', shell=True)
    if os.path.getsize(seq_file_dire+'all_Infor.txt') != 0:
        Information = pd.read_table(seq_file_dire+'all_Infor.txt', header=None)
        Information.columns = ['name','start','end','length','att','att_name','score','condition','dire','VFDB','SARG','Abi']
        Information['start'] = Information[['start']].astype(int)
        Information['end'] = Information[['end']].astype(int)
        Information = Information.sort_values(by=['name','start','end'], axis=0)
        Information = Information.drop_duplicates(inplace=False)
        for i in range(len(Information)):
            for j in range(i+1,len(Information)):
                if Information.iloc[i,0] != Information.iloc[j,0]:
                    break
                if Information.iloc[i,8] == Information.iloc[j,8] and Information.iloc[i,7] != 'uncomplete-int' and Information.iloc[j,7] != 'uncomplete-int':
                    if (int(Information.iloc[i,8]) == -1 and int(Information.iloc[j,1]) != int(Information.iloc[i,1]) and int(Information.iloc[j,1])-int(Information.iloc[i,1])<2500) or (int(Information.iloc[i,8]) == 1 and int(Information.iloc[j,1]) == int(Information.iloc[i,1])):
                        Information.iloc[j,4] = 'na'
                    if (int(Information.iloc[i,8]) == 1 and int(Information.iloc[j,1]) != int(Information.iloc[i,1]) and int(Information.iloc[j,2])-int(Information.iloc[i,2])<2500) or (int(Information.iloc[i,8]) == -1 and int(Information.iloc[j,1]) == int(Information.iloc[i,1])):
                        Information.iloc[i,4] = 'na'
        os.remove(seq_file_dire+'all_Infor.txt')
        Information = Information.drop(['dire'],axis = 1)
        Information.to_csv(seq_file_dire+'all_Infor.txt', sep='\t',index=False)
    time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]
    subprocess.run('echo "\033[33m ['+time+'] \033[0m" "\033[32mThank you for using PICI finder\033[0m"', shell=True)

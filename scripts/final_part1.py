import os
import time
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

def orf_hmm_prokka_diamond_find(filename,seq_file_name,threads,filepath,PICI_finder_path):
    if os.path.exists(seq_file_name + '_prokka'):
        return 0
    time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]
    cpu = int(seq_file_name.split('_')[-1])+1
    if int(threads) == 1:
        subprocess.run('echo "\033[33m ['+time+'] \033[0m" prokka finding ......', shell=True)
    if int(threads) > 1:
        subprocess.run('echo "\033[33m ['+time+'] \033[0m" progress '+str(cpu)+' prokka finding ......', shell=True)
    subprocess.run('prokka ' + filename + ' --outdir ' + seq_file_name + '_prokka' + ' --prefix PICI --centre X --compliant --metagenome --kingdom Bacteria --quiet', shell=True)
    time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]
    if int(threads) == 1:
        subprocess.run('echo "\033[33m ['+time+'] \033[0m" prokka finish', shell=True)
    if int(threads) > 1:
        subprocess.run('echo "\033[33m ['+time+'] \033[0m" progress '+str(cpu)+' prokka finish', shell=True)
    time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]
    if int(threads) == 1:
        subprocess.run('echo "\033[33m ['+time+'] \033[0m" prodigal finding ......', shell=True)
    if int(threads) > 1:
        subprocess.run('echo "\033[33m ['+time+'] \033[0m" progress '+str(cpu)+' prodigal finding ......', shell=True)
    subprocess.run('prodigal -a ' + seq_file_name + '_orf.txt' + ' -g 11 -p meta -q -i ' + filename + ' -o out.txt', shell=True)
    time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]
    if int(threads) == 1:
        subprocess.run('echo "\033[33m ['+time+'] \033[0m" prodigal finish', shell=True)
    if int(threads) > 1:
        subprocess.run('echo "\033[33m ['+time+'] \033[0m" progress '+str(cpu)+' prodigal finish', shell=True)
    if os.path.getsize(seq_file_name + '_orf.txt') == 0:
        return 0
    time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]
    if int(threads) == 1:
        subprocess.run('echo "\033[33m ['+time+'] \033[0m" hmmsearch finding ......', shell=True)
    if int(threads) > 1:
        subprocess.run('echo "\033[33m ['+time+'] \033[0m" progress '+str(cpu)+' hmmsearch finding ......', shell=True)
    p_list = []
    hmm1 = subprocess.Popen('hmmsearch --cpu 32 -o out.txt --domtblout ' + seq_file_name + '_hmmsearch.txt --noali --domE 1e-5 '+PICI_finder_path+'/databases/new_all_total.hmm ' + seq_file_name + '_orf.txt', shell=True)
    p_list.append(hmm1)
    hmm2 = subprocess.Popen('hmmsearch --cpu 32 -o out.txt --domtblout ' + seq_file_name + '_T4SS_hmm.txt --noali --domE 1e-5 '+PICI_finder_path+'/databases/T4SS_2.hmm ' + seq_file_name + '_orf.txt', shell=True)
    p_list.append(hmm2)
    hmm3 = subprocess.Popen('hmmsearch --cpu 32 -o out.txt --domtblout ' + seq_file_name + '_bac_hmm.txt --noali --domE 1e-5 '+PICI_finder_path+'/databases/Pfam-A-Bacteria.hmm ' + seq_file_name + '_orf.txt', shell=True)
    p_list.append(hmm3)
    for item in p_list:
        item.wait()
    subprocess.run('cat ' + seq_file_name + '_hmmsearch.txt ' + seq_file_name + '_T4SS_hmm.txt > ' + seq_file_name + '_all_hmm.txt', shell=True)
    time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]
    if int(threads) == 1:
        subprocess.run('echo "\033[33m ['+time+'] \033[0m" hmmsearch finish', shell=True)
    if int(threads) > 1:
        subprocess.run('echo "\033[33m ['+time+'] \033[0m" progress '+str(cpu)+' hmmsearch finish', shell=True)
    n = 0
    i = 0
    for i in range(1, len(open(seq_file_name + '_all_hmm.txt').readlines())):
        if linecache.getline(seq_file_name + '_all_hmm.txt', i).startswith('#'):
            n += 1
    if n == i:
        return 0
    time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]
    if int(threads) == 1:
        subprocess.run('echo "\033[33m ['+time+'] \033[0m" diamond finding ......', shell=True)
    if int(threads) > 1:
        subprocess.run('echo "\033[33m ['+time+'] \033[0m" progress '+str(cpu)+' diamond finding ......', shell=True)
    subprocess.run('diamond blastp --db '+PICI_finder_path+'/databases/all_diamond.dmnd --query ' + seq_file_name + '_orf.txt' + ' --out ' + seq_file_name + '_diamond_blast.txt --outfmt 6 --sensitive -k 1 --id 30 --quiet --threads 4', shell=True)
    subprocess.run('diamond blastp --db '+PICI_finder_path+'/databases/all_orf.dmnd --query ' + seq_file_name + '_orf.txt' + ' --out ' + seq_file_name + '_phage_diamond_blast.txt --outfmt 6 --sensitive -k 1 --id 30 --quiet --threads 4',shell=True)
    time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]
    if int(threads) == 1:
        subprocess.run('echo "\033[33m ['+time+'] \033[0m" diamond finish', shell=True)
    if int(threads) > 1:
        subprocess.run('echo "\033[33m ['+time+'] \033[0m" progress '+str(cpu)+' diamond finish', shell=True)
    
def orf_file_deal(seq_file_name,fasta_name):#将orf的数据变成数据框
    orf_file = pd.read_table(seq_file_name + '_orf.txt', header=None)
    orf_name_list = []
    diamond_list = []
    name_list = []
    gene_ano = ''
    flag = 0
    with open(seq_file_name + '_orf.txt') as orf_file_outcome:
        for line in orf_file_outcome.readlines():
            if line.startswith('>'):
                line = line.split('#')  # 以空格为分隔符
                orf_name_list.append(line)
    orf_name_data = pd.DataFrame(orf_name_list)
    orf_name_data[4] = orf_name_data[4].str.split(';',1).str[0]
    orf_name_data[4] = orf_name_data[4].str.split('ID=',1).str[1]
    orf_name_data[5] = orf_name_data[4].str.split('_',1).str[0]
    orf_name_data[6] = orf_name_data[4].str.split('_',1).str[1]
    orf_name_data[0] = orf_name_data[0].str.split('>',1).str[1]
    orf_name_data[7] = orf_name_data[0].str.rstrip()
    orf_name_data[0] = orf_name_data[0].str.rsplit('_',1).str[0]
    orf_name_data.columns = ['name','start', 'end', 'dir','mix', 'seq_number', 'gene_number','all_name']
    if os.path.getsize(seq_file_name + '_diamond_blast.txt') != 0:
       with open(seq_file_name + '_diamond_blast.txt') as diamond_outcome:
           for line in diamond_outcome.readlines():
               gene_ano = line.split()[1]
               if ((line.split()[1] != 'phage_resistance' or (line.split()[1] == 'phage_resistance' and int(line.split()[3]) < 300)) and line.split()[1] != 'stl_gene' and line.split()[1] != 'PICI_gene' and line.split()[1] != 'tail_gene' and line.split()[1] != 'head_gene' and line.split()[1] != 'lysis_gene' and line.split()[1] != 'capsid_gene' and float(line.split()[2]) < 90):
                   gene_ano = 'na'
               if line.split()[1] == 'tail_gene' and float(line.split()[2]) < 40:
                   gene_ano = 'low_tail_gene'
               name_list.append(gene_ano)
               diamond_list.append(line.split())
    diamond_list = pd.DataFrame(diamond_list)
    diamond_list[0] = diamond_list[0].str.rstrip()
    name_list = pd.DataFrame(name_list)
    all_diamond = pd.concat([diamond_list.iloc[:, [0]], name_list.iloc[:, [0]]], axis=1)
    all_diamond.columns = ['all_name','gene_name']
    orf = pd.merge(orf_name_data,all_diamond,on='all_name',how='outer')
    orf_table = pd.concat([orf.iloc[:, [0]], orf.iloc[:, [6]]], axis=1)
    orf_table = pd.concat([orf_table, orf.iloc[:, [1]]], axis=1)
    orf_table = pd.concat([orf_table, orf.iloc[:, [2]]], axis=1)
    orf_table = pd.concat([orf_table, orf.iloc[:, [3]]], axis=1)
    orf_table = pd.concat([orf_table, orf.iloc[:, [8]]], axis=1)
    orf_table = pd.concat([orf_table, orf.iloc[:, [5]]], axis=1)
    orf_table.fillna('na',inplace = True)
    return orf_table

def hmm_file_deal(seq_file_name,fasta_name):#将hmm的数据取最优E值，并且生存数据框
    hmm_list = []
    hmm_bac_list = []
    hmm_site = []
    with open(seq_file_name + '_all_hmm.txt') as hmmsearch_outcome:
        for line in hmmsearch_outcome.readlines():
            if line.startswith('#'):
                continue
            else:
                line = line.split()  # 以空格为分隔符
                hmm_list.append(line)
    with open(seq_file_name + '_bac_hmm.txt') as hmmsearch_bac_outcome:
        for line in hmmsearch_bac_outcome.readlines():
            if line.startswith('#'):
                continue
            else:
                line = line.split()  # 以空格为分隔符
                hmm_bac_list.append(line)
    hmm_list = pd.DataFrame(hmm_list)
    hmm_bac_list = pd.DataFrame(hmm_bac_list)
    hmm_bac_list.to_csv(seq_file_name + '_bac_hmmsearch.txt', sep='\t', index=False)
    hmm_list[2] = hmm_list[0].str.rsplit('_',1).str[1]
    hmm_list[0] = hmm_list[0].str.rsplit('_',1).str[0]
    hmm_list[5] = hmm_list[29].str.split(';',1).str[0]
    hmm_list[5] = hmm_list[5].str.split('ID=',1).str[1]
    hmm_list[5] = hmm_list[5].str.split('_',1).str[0]
    hmm_list[[5]] = hmm_list[[5]].astype(int)
    hmm_list[[2]] = hmm_list[[2]].astype(int)
    hmm_list[[11]] = hmm_list[[11]].astype(float)
    hmm_list = hmm_list.sort_values(by=[5,2,11], axis=0)
    hmm_list['temp'] = '_'
    hmm_list['new_name'] = hmm_list[11].map(str) + hmm_list['temp']
    hmm_list['new_name'] = hmm_list['new_name'].map(str) + hmm_list[12]
    hmm_list['new_name'] = hmm_list['new_name'].map(str) + hmm_list['temp']
    hmm_list['new_name'] = hmm_list['new_name'].map(str)  +hmm_list[6]
    hmm_list[11] = hmm_list['new_name']
    hmm_list.to_csv(seq_file_name + '_hmm_list', sep='\t', header = None ,index=False)
    list = [0, 11, 23, 25, 27, 2, 5, 19, 20]
    hmm_file_deal = pd.concat([hmm_list.iloc[:, [3]], hmm_list.iloc[:, [4]]], axis=1)
    for i in list:
        hmm_file_deal = pd.concat([hmm_file_deal, hmm_list.iloc[:, [i]]], axis=1)
    hmm_file_deal.dropna(axis=0, how='any', inplace=True)
    hmm_file_deal.columns = ['structure_name', 'structure_id', 'seq_name', 'E_value', 'start', 'end', 'dire', 'site','seq_number','domain_start','domain_end']
    return hmm_file_deal

def prokka_file_deal(seq_file_name,fasta_name):
    prokka_ext = []
    prokka_gene_ann = []
    with open(seq_file_name + '_prokka/' + 'PICI.gff') as prokka:
        for line in prokka.readlines():
            if line.startswith('gnl|X|') and line.split()[1] != 'prokka' and line.split()[2] != 'gene':
                prokka_ext.append(line.split('\t'))
    prokka_ext = pd.DataFrame(prokka_ext)
    prokka_ext[8] = prokka_ext[8].str.rstrip('\n')
    prokka_ext[8].to_csv(seq_file_name + '_prokka_ext_8', sep='\t',header = None, index=False)
    n = -1
    with open(seq_file_name + '_prokka_ext_8') as prokka_file:
        for line in prokka_file.readlines():
            n = n+1
            if 'note=CRISPR with' in line:
                prokka_ext.iloc[n,2] = 100
                continue
            if 'Name' in line:
                name_site=line.index('Name')
                name_str = line[:name_site]
                name_number = name_str.count(';')
                gene_name = line.split(';',name_number+1)[-2]
                prokka_ext.iloc[n,1] = gene_name.split('Name=')[1]
            else:
                prokka_ext.iloc[n,1] = 'na'
            produce_site=line.index('product')
            produce_str = line[:produce_site]
            produce_number = produce_str.count(';')
            gene_produce = line.split(';',produce_number+1)[-2]
            if 'locus_tag' in gene_produce:
                gene_produce = line.split(';',produce_number+1)[-1]
            if 'note=' in gene_produce:
                gene_produce = line.split(';',produce_number+1)[-1]
            prokka_ext.iloc[n,2] = gene_produce.split('product=')[1].rstrip('\n')
    prokka_ext = prokka_ext.replace('+', '1')
    prokka_ext = prokka_ext.replace('-', '-1')
    prokka_ext[7] = prokka_ext[0].str.split('_',1).str[1]
    prokka_ext = prokka_ext[~prokka_ext[2].isin([100])]
    list = [4, 6, 1, 2]
    prokka = pd.concat([prokka_ext.iloc[:, [7]], prokka_ext.iloc[:, [3]]], axis=1)
    for i in list:
        prokka = pd.concat([prokka, prokka_ext.iloc[:, [i]]], axis=1)
    prokka.columns = ['seq','start', 'end','dir', 'gene', 'protein']
    return prokka
    
if __name__ == "__main__":
    filename = sys.argv[1]
    seq_file_name = filename.rsplit(".",1)[0]
    cpu = int(seq_file_name.split('_')[-1])+1
    threads = sys.argv[2]
    threads_number = sys.argv[3]
    filepath = filename.rsplit('/', 1)[0]
    filepath = filepath.rsplit('/', 1)[0]
    PICI_finder_path = os.path.split(os.path.realpath(__file__))[0]
    PICI_finder_path = PICI_finder_path.rsplit('/', 1)[0]
    print(PICI_finder_path)
    n = 0
    name_list = []
    name_count = []
    name_length = []
    with open(filename) as fasta:
        for line in fasta.readlines():
            if line.startswith('>'):
                name_list.append(line.split('>')[1].rstrip('\n'))
                n += 1
                name_count.append(n)
            else:
                name_length.append(len(line))
    fasta_name = pd.concat([pd.DataFrame(name_list), pd.DataFrame(name_count)], axis=1)
    fasta_name = pd.concat([fasta_name, pd.DataFrame(name_length)], axis=1)
    fasta_name.to_csv(seq_file_name + '_fasta_name.txt', sep='\t', index=False)
    orf_hmm_prokka_diamond_find(filename,seq_file_name,threads,filepath,PICI_finder_path)
    time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]
    if int(threads) == 1:
        subprocess.run('echo "\033[33m ['+time+'] \033[0m" sequence information collection ......', shell=True)
    if int(threads) > 1:
        subprocess.run('echo "\033[33m ['+time+'] \033[0m" progress '+str(cpu)+' sequence information collection ......', shell=True)
    orf_table = orf_file_deal(seq_file_name,fasta_name)
    hmm_table = hmm_file_deal(seq_file_name,fasta_name)
    prokka_table = prokka_file_deal(seq_file_name,fasta_name)
    orf_table['temp'] = '_'
    orf_table['new_name'] = orf_table['seq_number'].map(str) + orf_table['temp']
    orf_table['new_name'] = orf_table['new_name'].map(str) + orf_table['start'].str.strip()
    orf_table['new_name'] = orf_table['new_name'].map(str) + orf_table['temp'].str.strip()
    orf_table['new_name'] = orf_table['new_name'].map(str) + orf_table['end'].str.strip()
    orf_table['new_name'] = orf_table['new_name'].map(str) + orf_table['temp'].str.strip()
    orf_table['new_name'] = orf_table['new_name'].map(str) + orf_table['dir'].str.strip()
    hmm_table['temp'] = '_'
    hmm_table['new_name'] = hmm_table['seq_number'].map(str) + hmm_table['temp']
    hmm_table['new_name'] = hmm_table['new_name'].map(str) + hmm_table['start'].str.strip()
    hmm_table['new_name'] = hmm_table['new_name'].map(str) + hmm_table['temp'].str.strip()
    hmm_table['new_name'] = hmm_table['new_name'].map(str) + hmm_table['end'].str.strip()
    hmm_table['new_name'] = hmm_table['new_name'].map(str) + hmm_table['temp'].str.strip()
    hmm_table['new_name'] = hmm_table['new_name'].map(str) + hmm_table['dire'].str.strip()
    prokka_table['temp'] = '_'
    prokka_table['new_name'] = prokka_table['seq'].map(str) + prokka_table['temp']
    prokka_table['new_name'] = prokka_table['new_name'].map(str) + prokka_table['start'].str.strip()
    prokka_table['new_name'] = prokka_table['new_name'].map(str) + prokka_table['temp'].str.strip()
    prokka_table['new_name'] = prokka_table['new_name'].map(str) + prokka_table['end'].str.strip()
    prokka_table['new_name'] = prokka_table['new_name'].map(str) + prokka_table['temp'].str.strip()
    prokka_table['new_name'] = prokka_table['new_name'].map(str) + prokka_table['dir'].str.strip()
    orf_hmm = pd.merge(orf_table,hmm_table,on='new_name',how='outer')
    orf_hmm_prokka = pd.merge(orf_hmm,prokka_table,on='new_name',how='outer')
    bac_table = pd.read_table(seq_file_name + '_bac_hmmsearch.txt',sep='\t')
    phage_diamond = pd.read_table(seq_file_name + '_phage_diamond_blast.txt', sep='\t', header=None)
    bac_table['bac_name'] = 'bac'
    new_bac_table = pd.concat([bac_table.iloc[:, [0]], bac_table.iloc[:, [30]]], axis=1)
    new_bac_table.columns = ['bac_phage_gene','bac_name']
    new_bac_table = new_bac_table.drop_duplicates(['bac_phage_gene'])
    phage_diamond[1] = 'phage'
    phage_diamond = phage_diamond.drop([2, 3, 4, 5, 6, 7, 8, 9, 10, 11], axis=1)
    phage_diamond.columns = ['bac_phage_gene','phage_name']
    orf_hmm_prokka['bac_phage_gene'] = orf_hmm_prokka['name'].map(str) + orf_hmm_prokka['temp_x'].str.strip()
    orf_hmm_prokka['bac_phage_gene'] = orf_hmm_prokka['bac_phage_gene'].map(str) + orf_hmm_prokka['gene_number'].str.strip()
    orf_hmm_prokka = pd.merge(orf_hmm_prokka,new_bac_table,on='bac_phage_gene',how='left')
    orf_hmm_prokka = pd.merge(orf_hmm_prokka,phage_diamond,on='bac_phage_gene',how='left')
    orf_hmm_prokka.fillna('na',inplace = True)
    orf_hmm_prokka.to_csv(seq_file_name + '_orf_hmm_prokka', sep='\t',index=False)
    orf_hmm_prokka.fillna('na',inplace = True)
    orf_hmm_prokka['seq_number_x'] = orf_hmm_prokka['new_name'].str.split('_').str[0]
    orf_hmm_prokka['start_x'] = orf_hmm_prokka['new_name'].str.split('_').str[1]
    orf_hmm_prokka['end_x'] = orf_hmm_prokka['new_name'].str.split('_').str[2]
    orf_hmm_prokka['dir_x'] = orf_hmm_prokka['new_name'].str.split('_').str[3]
    list = [3, 1, 4, 6, 7, 9, 10, 12, 25, 26, 5, 18, 19,29,30]
    orf_hmm_file_deal = pd.concat([orf_hmm_prokka.iloc[:, [0]], orf_hmm_prokka.iloc[:, [2]]], axis=1)
    for i in list:
        orf_hmm_file_deal = pd.concat([orf_hmm_file_deal, orf_hmm_prokka.iloc[:, [i]]], axis=1)
    orf_hmm_file_deal[['seq_number_x']] = orf_hmm_file_deal[['seq_number_x']].astype(int)
    orf_hmm_file_deal[['start_x']] = orf_hmm_file_deal[['start_x']].astype(int)
    orf_hmm_file_deal = orf_hmm_file_deal.sort_values(by=['seq_number_x','start_x','E_value'], axis=0)
    orf_hmm_file_deal.to_csv(seq_file_name + '_total_table.txt',sep='\t',index=False)
    os.remove(seq_file_name + '_all_hmm.txt')
    os.remove(seq_file_name + '_bac_hmmsearch.txt')
    os.remove(seq_file_name + '_bac_hmm.txt')
    os.remove(seq_file_name + '_diamond_blast.txt')
    os.remove(seq_file_name + '_hmm_list')
    os.remove(seq_file_name + '_hmmsearch.txt')
    os.remove(seq_file_name + '_orf_hmm_prokka')
    os.remove(seq_file_name + '_phage_diamond_blast.txt')
    os.remove(seq_file_name + '_prokka_ext_8')
    os.remove(seq_file_name + '_T4SS_hmm.txt')
    shutil.rmtree(seq_file_name + '_prokka')
    time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]
    if int(threads) == 1:
        subprocess.run('echo "\033[33m ['+time+'] \033[0m" sequence information finish', shell=True)
    if int(threads) > 1:
        subprocess.run('echo "\033[33m ['+time+'] \033[0m" progress '+str(cpu)+' sequence information collecion finish', shell=True)
    subprocess.run('python3 '+PICI_finder_path+'/scripts/final_part3.py '+seq_file_name+' '+filepath+' '+str(threads), shell=True)


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


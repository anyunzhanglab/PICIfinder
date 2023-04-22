#! /usr/bin/env python3
# Author: Yuxuan Feng

# PICI finder v1.0.0
# PICI Identification From Genome
# Release date: 10 6 2022
# Usage: $ python3 PICI_finder.py --seqfile <input_file> [options]

############################### Imports  #######################################

import os
import math
import traceback
import subprocess
import linecache
import pandas as pd
import click
from datetime import datetime
from Bio.SeqIO.FastaIO import SimpleFastaParser

############################### Options  #######################################
@click.command()
@click.option('-i','--seqfile', type=click.Path(resolve_path=True),required=True, help='input fasta file')
@click.option('-t','--cpu',type=int,default=1, help='number of parallel PICI_run runs]')
@click.option('-l','--length',type=int,default=3000, help='length in basepairs to limit input sequences [default=1000, can increase but not decrease]')
@click.option('-o','--out',type=click.Path(dir_okay=True,writable=True,resolve_path=True), default='.',help="path to deposit output folder and temporary files, will create if doesn't exist [default= working directory]")

def cpu_split(seqfile,cpu,length,out):
    out = out+'/'
    PICI_finder_path = os.path.split(os.path.realpath(__file__))[0]

    if not os.path.exists(out):
        subprocess.run('mkdir ' + out, shell=True)
    
    try:
        filename = str(seqfile.rsplit('/')[-1])
        filename_resuf = str(filename.rsplit(".",1)[0])
        filepath = str(seqfile.rsplit('/')[0])
    except Exception:
        filename = str(seqfile)

    sequences = 0
    with open(str(seqfile),'r') as infile:
        file_seq_name = open(out + str(filename_resuf) + '_seq_name.txt', 'w')
        complete_file = open(out + str(filename_resuf)+'_sequence.txt','w')
        for name, seq in SimpleFastaParser(infile):
            if len(str(seq)) >= int(length):
                list = [' ', '.', '/', ',', ')', '(', ':', '|', ';', '-', ';', '\'', '#', '=', '*']
                for line in list:
                    name = name.replace(line, '_')
                name = '>'+name
                file_seq_name.write(str(name) + '\n')
                complete_file.write(str(name) + '\n')
                complete_file.write(str(seq) + '\n')
                sequences += 1
        file_seq_name.close()
        complete_file.close()
        
    try:
        lines = math.ceil(int(sequences)/cpu)
        parallel = 'split -l ' + str(lines) + ' ' + out + str(filename_resuf)+'_seq_name.txt' + ' ' + out + str(filename_resuf)+'.parallel-runs_'
        subprocess.run(parallel, shell=True)
    except Exception as e:
        print(e.__class__.__name__)
        print(e)
        print('It may be no sequence greater than '+str(length)+' in length in the file.')
        subprocess.run('rm ' + out + str(filename_resuf)+'_seq_name.txt', shell=True)
        subprocess.run('rm ' + out + str(filename_resuf)+'.txt', shell=True)
        subprocess.run('rm ' + out + str(filename_resuf)+'.parallel-runs_', shell=True)
        return 0

    if not os.path.exists(PICI_finder_path+'/databases/all_diamond.dmnd'):
        subprocess.run('diamond makedb --in '+PICI_finder_path+'/databases/all_diamond2.fasta --db '+PICI_finder_path+'/databases/all_diamond.dmnd', shell=True)
    if not os.path.exists(PICI_finder_path+'/databases/all_orf.dmnd'):
        subprocess.run('diamond makedb --in '+PICI_finder_path+'/databases/all_orf.ffa --db '+PICI_finder_path+'/databases/all_orf.dmnd', shell=True)
    
    parallels = 'ls ' + out + str(filename_resuf) + '.parallel-runs_*'
    parallels_shell = subprocess.check_output(parallels, shell=True)
    parallels1 = str(parallels_shell.strip()).split("'")[1]
    parallels2 = parallels1.split("\\n")

    n = 0
    p_list = []
    while n < len(parallels2):
        fna = open(str(parallels2[n]).rsplit(".",1)[0]+'_'+str(n)+'.ffa','a')
        for i in range(1,len(open(parallels2[n]).readlines())+1):
            for j in range(1,len(open(out + str(filename_resuf)+'_sequence.txt').readlines())+1):
                if linecache.getline(parallels2[n],i) == linecache.getline(out + str(filename_resuf)+'_sequence.txt',j):
                    fna.write(linecache.getline(out + str(filename_resuf)+'_sequence.txt',j))
                    fna.write(linecache.getline(out + str(filename_resuf)+'_sequence.txt', j+1))
                    continue
        fna.close()
        subprocess.run('rm ' + parallels2[n], shell=True)
        n += 1

    parallels = 'ls ' + out + str(filename_resuf) + '*.ffa'
    parallels_shell = subprocess.check_output(parallels, shell=True)
    parallels1 = str(parallels_shell.strip()).split("'")[1]
    parallels2 = parallels1.split("\\n")
    
    os.remove(out + str(filename_resuf)+'_seq_name.txt')
    os.remove(out + str(filename_resuf)+'_sequence.txt')
    n = 0
    p_list = []
    while n < len(parallels2):
        judge = 'python3 '+PICI_finder_path+'/scripts/final_part1.py ' + parallels2[n] + ' ' + str(cpu) + ' ' + str(n) 
        level1 = subprocess.Popen(judge, shell=True)
        p_list.append(level1)
        n += 1
    
    for item in p_list:
        item.wait()

    subprocess.run('python3 '+PICI_finder_path+'/scripts/final_part2.py '+out, shell=True)

if __name__ == '__main__':
    time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]
    subprocess.run('echo "\033[33m ['+time+'] \033[0m" "\033[31mPICI finder start !!!\033[0m"', shell=True)
    cpu_split()

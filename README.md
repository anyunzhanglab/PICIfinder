# PICIfinder
PICIfinder is a tool that provide utilities to detect PICI directly from the Gram-positive bacteria genome.It is mainly based on the conserved pattern features of PICI, and uses hmmsearch to compare the input genome sequences to find the fragments that match the PICI pattern.The detection algorithm consists of four parts: (1) Prediction of candidate fragments based on co-localization of gene modules.(2) Distinguishing false positives of PICI candidate fragments.(3).Range determination of PICI by random forest model.(4) Determination of the grade of PICI.
Quick Start
Download and install

git clone https://github.com/anyunzhanglab/PICIfinder.git
cd PICIfinder
.....

Run PICIfinder 

PICIfinder.py \
	-i my.fasta \
	-l 1000 \
	-t 10 \
 	-o /OUTPATH \
Installation
Requirements
Before running PICIfinder, you need to install python packages and some tools.
PICIfinder is written in Python and depends on the following packages:Bio.SeqIO.FastaIO，sklearn.ensemble，sklearn.model_selection
	Python packages
	Bio.SeqIO
	FastaIO
	sklearn.ensemble
	sklearn.model_selection
	.Tools
	prokka
	prodigal
	hmmsearch
	diamond
PICIfinder Usage
Usage: PICI_finder.py [OPTIONS]

Options:
  -i, --seqfile PATH    input fasta file  [required]
  -t, --cpu INTEGER     number of parallel PICI_run runs]
  -l, --length INTEGER  length in basepairs to limit input sequences
                        [default=1000, can increase but not decrease]
  -o, --out PATH        path to deposit output folder and temporary files,
                        will create if doesn't exist [default= working
                        directory]
  --help                Show this message and exit.
Example
	PICIfiner.py -i example/SaPI1.fasta -o outpath/SaPI1

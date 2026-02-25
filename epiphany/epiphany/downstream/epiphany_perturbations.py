# 1. Load packages
import torch
import torch.nn as nn
import torch.nn.functional as F
import torchvision.transforms as transforms
from torch import randn
from torch.nn import MSELoss
import torch.optim as optim
from torch.optim import Adam
from torch.utils.data import DataLoader
import numpy as np
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt
import random
import pandas as pd
import seaborn as sns
import time
import pickle
from datetime import datetime
from torch.autograd import Variable
import gzip
import sys
import os 
from sklearn.decomposition import TruncatedSVD, PCA
import pyBigWig
import argparse

from utils.model_architecture_util import Net_wider,restore
from utils.generate_perturbations_util import data_load, prediction
# Create ArgumentParser object
parser = argparse.ArgumentParser(description="Parse input arguments for cell type.")
parser.add_argument('--c', type=str, required=False, default="GM12878", help="Cell type")
parser.add_argument('--path', type=str, required=True, help="output path")
# Parse the arguments
args = parser.parse_args()
cell_type = args.c
output_path= args.path

# 2. Load data - part 2
chrom_list = ["chr"+str(i) for i in range(1,23)] #for human hg38
length_list = [240000000,240000000,190000000,190000000,180000000,170000000,150000000,
            140000000,130000000,130000000,130000000,130000000,110000000,100000000,100000000,90000000,
            80000000, 80000000,50000000,60000000,40000000,50000000]
chrom_len_dict = dict(zip(chrom_list,length_list))

#Load model 
wsize = 14000
net = Net_wider(input_channels=5,window_size=wsize)
model_path = '/scratch/gent/vo/000/gvo00027/projects/CBIGR/24NHT_master_thesis/code/epiphany/models/pretrained_GM12878_5kb.pt_model'
restore(net,model_path)
net.eval()

for chrom, length in chrom_len_dict.items():
    # Extract epigenomic data from downloaded tracks
    bwfile_dir = f"/scratch/gent/vo/000/gvo00027/projects/CBIGR/24NHT_master_thesis/code/epiphany/bigWig/{cell_type}"
    input_tracks = data_load(chrom=chrom, bwfile_dir=bwfile_dir, distances=100, cell_type=cell_type)
    start_offset=1000000
    min_slice= 2000000
    effective_length = length - start_offset
    n=10
    slice_size = effective_length // n    
    if slice_size < min_slice:
        slice_size = min_slice
    else:
        slice_size = (slice_size//1000000)*1000000
    start = start_offset
    while start < (length-2000000):
        bin0 = start+1000000 - 2500
        bin1 = bin0 + 5000
        #bin0=1997500
        #bin1=2002500
        prediction(chrom,net,cell_type,input_tracks,bin0,bin1, save_txt=f"{output_path}/temp_{chrom}_{start}.txt.gz")
        prediction(chrom,net,cell_type,input_tracks,bin0,bin1,perturb_type="deletion", save_txt=f"{output_path}/temp_del_{chrom}_{start}.txt.gz")
        start+=slice_size
        print(f"Processed {cell_type} {chrom} {start}")
        

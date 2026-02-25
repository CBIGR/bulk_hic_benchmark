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
from utils.generate_predictions_util import *
import argparse
from utils.model_architecture_util import * 

parser = argparse.ArgumentParser(description="Parse input arguments for cell type.")
parser.add_argument('--model_path', type=str, required=True,  help="model_path")
parser.add_argument('--gt_file', type=str, required=True, help="normalized groundtruth path")

parser.add_argument('--bwfile_dir', type=str, required=True,  help="directory of bw files")
parser.add_argument('--dataset', type=str, required=True, help="Cell type")
parser.add_argument('--output_path', type=str, required=True, help="Path to save file")
parser.add_argument('--chr', type=str, required=True, help="Path to save file")

# Parse the arguments
args = parser.parse_args()


model_path = args.model_path
chrom = args.chr
gt_file = args.gt_file
bwfile_dir = args.bwfile_dir
dataset = args.dataset
output_path = args.output_path
chr=args.chr

#Load model 
wsize = 14000
net = Net(input_channels=5,window_size=wsize)
restore(net,model_path)
net.eval()

chrom = "chr"+chr
print(chrom,datetime.now())
results_generation(chrom = chrom, net=net, 
                    cell_type = dataset, 
                    bwfile_dir = bwfile_dir,
                    submatrix_location = f"{output_path}/{chrom}_intermediate_matrices.txt", assemble_matrix_location = f"{output_path}/{chrom}_assembled_chromosome.txt",
                    ground_truth_file = gt_file, ground_truth_location = f"{output_path}/{chrom}_ground_truth_corresponding_location.txt",window_size=14000,
                    seq_length = 200,resolution_hic = 10000) #normcounts, zvalue, zfull
     
     
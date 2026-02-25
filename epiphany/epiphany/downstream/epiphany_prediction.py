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
torch.set_default_tensor_type(torch.DoubleTensor)
import pyBigWig
from utils.model_architecture_util import Net, restore
from utils.generate_predictions_util import results_generation

chrom_list = ["chr"+str(i) for i in range(1,23)] #for human hg38
length_list = [248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,
               138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,
               83257441,80373285,58617616,64444167,46709983,50818468]
chrom_len_dict = dict(zip(chrom_list,length_list))

#Load model 
wsize = 14000
net = Net(input_channels=5,window_size=wsize)
model_path = '/scratch/gent/vo/000/gvo00027/projects/CBIGR/24NHT_master_thesis/code/epiphany/models/pretrained_GM12878.pt_model'
restore(net,'/content/epiphany/pretrained/pretrained_GM12878.pt_model')
net.eval()

#Generate predictions for a single chromosome
chrom = "chr3"
bwfile_dir = "/scratch/gent/vo/000/gvo00027/projects/CBIGR/24NHT_master_thesis/code/epiphany/bigWig/GM12878"
ground_truth_file = "bwfile_dir"
gt_path = "/scratch/gent/vo/000/gvo00027/projects/CBIGR/24NHT_master_thesis/code/epiphany/ground_truth/Apr8_GM12878_chr3_normcounts_diagonal.txt"
print(chrom,datetime.now())
results_generation(chrom = chrom, net=net, 
                    cell_type = "GM12878", 
                    bwfile_dir = bwfile_dir,
                    submatrix_location = "./intermediate_matrices.txt", assemble_matrix_location = "./assembled_chromosome.txt",
                    ground_truth_file = gt_path, ground_truth_location = "/content/ground_truth_corresponding_location.txt", 
                    window_size = wsize) #normcounts, zvalue, zfull

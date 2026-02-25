import numpy as np
import pandas
import argparse
import pandas as pd

# Initialize parser
parser = argparse.ArgumentParser(description="A simple script with argparse")

# Add arguments
parser.add_argument("--path", help="Dir to the prediction file")
parser.add_argument("--dataset",help="Dataset")
parser.add_argument("--bedpe_path_save", help="Dir to save bedpe file for generate HIC")

# Parse arguments
args = parser.parse_args()
CHROMLENGTH = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
}

res=10000
dataset = args.dataset
path = args.path
bedpe_path_save = args.bedpe_path_save
bedpes = []
window=2000000
b=200
for c in range(1,23):
    start_mat = 1000000
    chr=f"chr{c}"
    file_pred=f"{path}/prediction_{dataset}_{chr}.npz"
    with open(file_pred,"rb") as f:
        mat = np.load(f)['arr_0']
    # convert back to df start, end, count
    starts = []
    ends = []
    norm_count = []
    for r in range(mat.shape[0]):
        starts += list(range(start_mat,(start_mat+window),res)) + [start_mat+window]*b
        ends += [start_mat+window]*b + list(range(start_mat+window+res,(start_mat+res+window*2),res))
        norm_count += mat[r].tolist() 
        start_mat += res 
    df_hic = pd.DataFrame(list(zip(starts, ends, norm_count)), columns=["start","end","norm_count"])
    df_hic['start'] = df_hic['start'].astype(int)
    df_hic['end'] = df_hic['end'].astype(int)
    df_hic = df_hic.groupby(['start', 'end'], as_index=False)['norm_count'].mean()
    df_hic.sort_values(by=['start','end'],inplace=True)
    bedpe = pd.DataFrame({
    "col0": [0] * len(df_hic),
    "chrom1": [c]*len(df_hic),
    "start1": df_hic['start'].values,
    "col01": [0] * len(df_hic),
    "col02": [0] * len(df_hic),
    "chrom2": [c]*len(df_hic),
    "start2": df_hic['end'].values,
    "col1": [1] * len(df_hic),
    "count": df_hic['norm_count'].values})

    bedpes += [bedpe]

merge_bedpe = pd.concat(bedpes, ignore_index=True)
merge_bedpe.to_csv(f"{bedpe_path_save}/{dataset}_{res}_norm.bedpe", sep="\t", header=False, index=False)

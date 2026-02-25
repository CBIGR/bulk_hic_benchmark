import os
import numpy as np
from scipy import sparse
from tqdm import tqdm 
import argparse
import pandas as pd
import cooler

parser = argparse.ArgumentParser(description='Convert npy to cool')
parser.add_argument('--dataset', required=True, help='dataset')
parser.add_argument('--datadir', required=True, help='directory contain data for convertion')
parser.add_argument('--output', required=True, help='directory to save bedpe file')
parser.add_argument('--resolution', required=True, help='resolution')
parser.add_argument('--window_bp', required=True, help='window size')

args = parser.parse_args()
dataset=args.dataset
data_dir=args.datadir
output_dir=args.output


resolution = int(args.resolution) # 8192  # bin size
window_bp = int(args.window_bp) # 2097152  # 2Mb windows
bins_per_window = window_bp // resolution  # 256

# Chromosome sizes (hg38 – adjust as needed)
chrom_sizes = {
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


bedpes = []
bedpes_norm = []
for chrom, chrom_size in chrom_sizes.items():
    print(f"\nProcessing {chrom}...")

    num_bins = chrom_size // resolution + 1
    global_matrix = sparse.lil_matrix((num_bins, num_bins), dtype=np.float32)

    files = sorted([f for f in os.listdir(data_dir) if f.startswith(chrom + "_") and f.endswith(".npz")])

    for fname in tqdm(files, desc="Merging submatrices"):
        # Extract start bin from filename
        try:
            # Expected format: chr1_00000000.npy
            start_bp = int(fname.replace(f"{chrom}_", "").replace(".npz", ""))
        except ValueError:
            print(f"Skipping malformed file: {fname}")
            continue
        start_bin = start_bp // resolution
        print(f"Processing bp:{start_bp} bin:{start_bin}")
        end_bin = start_bin + bins_per_window

        # Skip if out of bounds
        if end_bin > num_bins:
            continue

        # Load submatrix
        submatrix = np.load(os.path.join(data_dir, fname))['data']

        # Insert into global matrix
        global_matrix[start_bin:end_bin, start_bin:end_bin] = submatrix

    # Convert to cooler-compatible sparse format
    global_matrix = global_matrix.tocoo()

    print(f"Done assembling global matrix {chrom}.")
    row, col, data = global_matrix.row, global_matrix.col, global_matrix.data
    count = np.exp(data)-1
    
    bins = pd.DataFrame({
    "chrom": [chrom[3:]] * num_bins,
    "start": [i * resolution for i in range(num_bins)],
    "end": [(i + 1) * resolution for i in range(num_bins)],
    })
    
    # Create BEDPE from bin positions
    bedpe = pd.DataFrame({
        "col0": [0] * len(row),
        "chrom1": bins.loc[row, "chrom"].values,
        "start1": bins.loc[row, "start"].values,
        "col01": [0] * len(row),
        "col02": [0] * len(row),
        "chrom2": bins.loc[col, "chrom"].values,
        "start2": bins.loc[col, "start"].values,
        "col1": [1] * len(row),
        "count": count,
    })
    bedpe.sort_values(by=['start1','start2'],inplace=True)
    
    # Create BEDPE from bin positions
    bedpe_norm = pd.DataFrame({
        "col0": [0] * len(row),
        "chrom1": bins.loc[row, "chrom"].values,
        "start1": bins.loc[row, "start"].values,
        "col01": [0] * len(row),
        "col02": [0] * len(row),
        "chrom2": bins.loc[col, "chrom"].values,
        "start2": bins.loc[col, "start"].values,
        "col1": [1] * len(row),
        "count": data,
    })
    
    bedpe_norm.sort_values(by=['start1','start2'],inplace=True)
    
    bedpes_norm += [bedpe_norm]
    bedpes += [bedpe]
    


    print(f"Done processing {chrom}.")

# Now that all chromosomes are processed, finalize the .mcool file
# Once all chromosomes have been added, save it as a .mcool file
print("Creating the final .mcool file with multiple chromosomes...")

# Optionally, you can add extra metadata, such as different resolutions or attributes, if necessary.

merge_bedpe = pd.concat(bedpes, ignore_index=True)
merge_bedpe.to_csv(f"{output_dir}/{dataset}_{resolution}_count.bedpe", sep="\t", header=False, index=False)

merge_bedpe = pd.concat(bedpes_norm, ignore_index=True)
merge_bedpe.to_csv(f"{output_dir}/{dataset}_{resolution}_norm.bedpe", sep="\t", header=False, index=False)
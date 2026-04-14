import numpy as np
import LASCA_pipeline
import cooler

import argparse

parser = argparse.ArgumentParser(description="Extract Hi-C matrix from .hic file")

parser.add_argument('--hic', required=True, help='Path to .hic file')
parser.add_argument('--res', type=int, default=10000, help='Resolution (bin size), e.g., 10000')
parser.add_argument('--distance',type=int, help='Maximum loop size')
parser.add_argument('--output', help='Directory to save output file')

args = parser.parse_args()

# Load Hi-C matrix
datapath=args.hic
resolution=args.res
distance_bins = args.distance
output=args.output

#import cool file
S=cooler.Cooler(f"{datapath}::resolutions/10000")

#list of chromosomes
chroms=S.chromnames

for c in chroms:
    
    #load HiC contacts
    #raw count
    mtx_raw_S=S.matrix(balance=False).fetch(c)
    #Balanced
    mtx_S=S.matrix(balance=True).fetch(c)
    
    # Define bin range
    start_bin = 0
    end_bin = mtx_S.shape[0] - 1
    
    # Run LASCA
    LASCA_pipeline.LASCA_processing(
        raw_mtx=mtx_raw_S,
        mtx=mtx_S,
        chr_name="chr"+str(c),
        output_bedpe_name=f"{output}/LASCA_{c}_{resolution}.bedpe",
        resolution=resolution,
        start_bin=start_bin,
        end_bin=end_bin,
        distance_bins=distance_bins,
        FDR=0.1,
        q=0.9,
        adjust_by_scale=False,
        min_cluster=3,
        filter_bins=2,
        q_value_trhd=0.05,
        Intensity=True,
        as_intervals=True,
        bin_coverage=0.25,
        save_qvalues_mtx=False,
        use_filter=False
    )
    
    del mtx_raw_S
    del mtx_S
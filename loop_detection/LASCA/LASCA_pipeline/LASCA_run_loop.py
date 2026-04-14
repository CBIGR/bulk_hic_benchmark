import numpy as np
import LASCA_pipeline
import hicstraw

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
hic = hicstraw.HiCFile(datapath)
start = 0


for c in range(1,23):
    
    end = hic.getChromosomes()[c].length
    
    # Load raw contacts (unbalanced)
    mzd_raw = hic.getMatrixZoomData(str(c), str(c), 'observed', 'NONE', 'BP', resolution)
    mtx_raw_S = mzd_raw.getRecordsAsMatrix(3000000, end, start, end)
    
    # Load balanced contacts (if available)
    mzd_norm = hic.getMatrixZoomData(str(c), str(c), 'observed', 'KR', 'BP', resolution)
    mtx_S = mzd_norm.getRecordsAsMatrix(3000000, end, start, end)
    
    # Define bin range
    start_bin = 0
    end_bin = mtx_S.shape[0] - 1
    
    # Run LASCA
    LASCA_pipeline.LASCA_processing(
        raw_mtx=mtx_raw_S,
        mtx=mtx_S,
        chr_name="chr"+str(c),
        output_bedpe_name=f"{output}/LASCA_chr{str(c)}_{resolution}.bedpe",
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
    del mzd_raw
    del mtx_raw_S
    del mzd_norm
    del mtx_S
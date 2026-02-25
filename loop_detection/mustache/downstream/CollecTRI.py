import pandas as pd
import decoupler as dc
import pyranges as pr
from pybedtools import BedTool
import pybedtools
import re
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description="Annotate Hi-C loops with enhancer and promoter information.")

# Input loop file (BEDPE)
parser.add_argument("--model", help="Model name")

# Enhancer annotation file (BED)
parser.add_argument("--path_save", help="Path to save")

# Promoter annotation file (BED)
parser.add_argument(
    "--path", help="Path to loop")
args = parser.parse_args()
loop_file=args.path
path_save = args.path_save
model = args.model

def has_tf_target(anchor1_genes, anchor2_genes, regulons):
    for gene1 in anchor1_genes:
        for gene2 in anchor2_genes:
            if ((regulons['source'] == gene1) & (regulons['target'] == gene2)).any() or \
               ((regulons['source'] == gene2) & (regulons['target'] == gene1)).any():
                return True
    return False

gene_bed = "gencode.v48.basic.annotation.gtf.gz"

loops = pd.read_csv(loop_file, sep="\t")
col_names=['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2', 'FDR', 'DETECTION_SCALE']
loops.columns = col_names
col_convert = ['chr1','x1','x2','chr2','y1', 'y2']
loops = loops[~loops['chr1'].isin(['X','Y'])]
loops[col_convert] = loops[col_convert].astype(int)
loops['chr1'] = 'chr' + loops['chr1'].astype('str')
loops['chr2'] = 'chr' + loops['chr2'].astype('str')

regulons = dc.get_collectri(organism='human', split_complexes=False)
regulons = regulons[['source', 'target']]

genes = pd.read_csv(gene_bed, sep="\t", header=None,  comment='#',compression='gzip')
col_names = ["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
genes.columns = col_names


# Filter to TSS/promoters of genes
promoters = genes[genes["feature"] == "transcript"]
promoters["TSS"] = promoters.apply(lambda x: x["start"] if x["strand"] == "+" else x["end"], axis=1)
promoters["start_promoter"] = promoters["TSS"] - 2000
promoters["end_promoter"] = promoters["TSS"] + 2000

# Clean up
promoter_bed = promoters[["chrom", "start_promoter", "end_promoter"]].copy()
promoter_bed.to_csv(f"{path_save}/promoters.bed", sep="\t", header=False, index=False)

genes = genes[genes['feature'] == 'gene']
genes['gene_id'] = genes['attribute'].apply(lambda x: x.split(';')[0][9:-1])
genes['gene_name'] = genes['attribute'].apply(lambda x: x.split(';')[2][12:-1])

anchor1 = loops[['chr1', 'x1', 'x2']].copy()
anchor1.columns = ['chrom', 'start', 'end']

anchor2 = loops[['chr2', 'y1', 'y2']].copy()
anchor2.columns = ['chrom', 'start', 'end']

anchor1_bed = BedTool.from_dataframe(anchor1)
anchor2_bed = BedTool.from_dataframe(anchor2)

genes_bed = BedTool.from_dataframe(genes[['chrom', 'start', 'end', 'gene_name']])

overlap1 = anchor1_bed.intersect(genes_bed, wa=True, wb=True).to_dataframe(names=[
    'chrom', 'start', 'end', 'gchrom', 'gstart', 'gend','gene_name'])
overlap2 = anchor2_bed.intersect(genes_bed, wa=True, wb=True).to_dataframe(names=[
    'chrom', 'start', 'end', 'gchrom', 'gstart', 'gend', 'gene_name'])

loop_genes = loops.copy()

loop_genes['anchor1_genes'] = overlap1.groupby([
    'chrom', 'start', 'end'])['gene_name'].apply(list).reindex(anchor1.set_index([
        'chrom', 'start', 'end']).index).reset_index(drop=True)
loop_genes['anchor2_genes'] = overlap2.groupby([
    'chrom', 'start', 'end'])['gene_name'].apply(list).reindex(anchor2.set_index([
        'chrom', 'start', 'end']).index).reset_index(drop=True)

loop_genes['anchor1_genes'] = loop_genes['anchor1_genes'].apply(lambda x: x if isinstance(x, list) else [])
loop_genes['anchor2_genes'] = loop_genes['anchor2_genes'].apply(lambda x: x if isinstance(x, list) else [])

target_to_tfs = defaultdict(set)
for _, row in regulons.iterrows():
    target_to_tfs[row['target']].add(row['source'])

def get_matched_tfs(genes_list):
    tfs = set()
    for gene in genes_list:
        tfs.update(target_to_tfs.get(gene, set()))
    return tfs

def matched_tfs_for_loop(row):
    tfs = get_matched_tfs(row['anchor1_genes']).union(get_matched_tfs(row['anchor2_genes']))
    return ';'.join(tfs) if tfs else None

loop_genes['matched_tfs'] = loop_genes.apply(matched_tfs_for_loop, axis=1)

loop_genes['is_matched_collecttri'] = loop_genes['matched_tfs'].apply(lambda x: 1 if x else 0)


# Optional: summarize
n_total = len(loop_genes)
n_with_tf_target = loop_genes['is_matched_collecttri'].sum()
print(f"Total loops: {n_total}")
print(f"Loops with TF-target interaction: {n_with_tf_target} ({(n_with_tf_target/n_total)*100:.2f}%)")

# ENHANCER
enhancer_bedpe = pd.read_csv("enhancer_K562.bed", sep='\t', header=None)
enhancer_bedpe.columns = ['enh_chr','enh_start','enh_end','value']


enh_gr = pr.PyRanges(pd.DataFrame({
    "Chromosome": enhancer_bedpe ['enh_chr'],
    "Start": enhancer_bedpe ['enh_start'],
    "End": enhancer_bedpe ['enh_end']
}))

anchor1 = pr.PyRanges(pd.DataFrame({
    "Chromosome": loop_genes['chr1'],
    "Start": loop_genes['x1'],
    "End": loop_genes['x2'],
    "loop_index": loop_genes.index
}))

anchor2 = pr.PyRanges(pd.DataFrame({
    "Chromosome": loop_genes['chr2'],
    "Start": loop_genes['y1'],
    "End": loop_genes['y2'],
    "loop_index": loop_genes.index
}))

anchor1_overlap = anchor1.join(enh_gr)
anchor2_overlap = anchor2.join(enh_gr)

# Get unique indices of overlapping anchors
anchor1_indices = set(anchor1_overlap.df['loop_index'])
anchor2_indices = set(anchor2_overlap.df['loop_index'])

loop_genes['anchor1_is_enhancer'] = loop_genes.index.isin(anchor1_indices)
loop_genes['anchor2_is_enhancer'] = loop_genes.index.isin(anchor2_indices)

# promoter
promoter_gr = pr.PyRanges(pd.DataFrame({
    "Chromosome": promoter_bed['chrom'],
    "Start": promoter_bed['start_promoter'],
    "End": promoter_bed['end_promoter']
}))


# Intersect
anchor1_overlap = anchor1.join(promoter_gr)
anchor2_overlap = anchor2.join(promoter_gr)

# Get unique indices of overlapping anchors
anchor1_indices = set(anchor1_overlap.df['loop_index'])
anchor2_indices = set(anchor2_overlap.df['loop_index'])

loop_genes['anchor1_is_promoter'] = loop_genes.index.isin(anchor1_indices)
loop_genes['anchor2_is_promoter'] = loop_genes.index.isin(anchor2_indices)

# chip atlas

# Assign column names (based on BED + extra columns)
tf_binding = pd.read_csv("K-562_TF_chipatlas.bed", sep= '\t', skiprows=1, header=None)
tf_binding.columns = ['chrom', 'start', 'end', 'attributes', 'score', 'strand', 'thickStart', 'thickEnd', 'rgb']

# Extract TF name from the "Name=" field in the attributes
def extract_tf(attr):
    match = re.search(r'Name=([^%;]+)', attr)
    return match.group(1).split('%')[0] if match else None

tf_binding['TF'] = tf_binding['attributes'].apply(extract_tf)

# Select useful columns for intersection
chipatlas_clean = tf_binding[['chrom', 'start', 'end', 'TF']]
chipatlas_clean.columns = ['Chromosome', 'Start', 'End', 'TF_name']
chipatlas_pr = pr.PyRanges(chipatlas_clean)
del tf_binding

# --- Intersect anchors with ChIP-Atlas ---
overlap1 = anchor1.join(chipatlas_pr)
overlap2 = anchor2.join(chipatlas_pr)

# --- Merge overlaps back with loop dataframe ---
# For anchor1
df1 = overlap1.df.groupby("loop_index")["TF_name"].apply(lambda x: ",".join(sorted(set(x)))).reset_index()
df1.rename(columns={"TF_name": "TF_anchor1"}, inplace=True)

# For anchor2
df2 = overlap2.df.groupby("loop_index")["TF_name"].apply(lambda x: ",".join(sorted(set(x)))).reset_index()
df2.rename(columns={"TF_name": "TF_anchor2"}, inplace=True)

# Merge back into loops DataFrame
loop_genes = loop_genes.merge(df1, left_index=True, right_on="loop_index", how="left").drop(columns="loop_index")
loop_genes = loop_genes.merge(df2, left_index=True, right_on="loop_index", how="left").drop(columns="loop_index")

# Optional: Add single column with all TFs
loop_genes["TFs_in_loop"] = loop_genes[["TF_anchor1", "TF_anchor2"]].fillna("").agg(",".join, axis=1)
loop_genes["TFs_in_loop"] = loop_genes["TFs_in_loop"].str.strip(",").str.replace(",,", ",").replace("", pd.NA)
loop_genes['is_matched_tf_chipatlas'] = loop_genes['TFs_in_loop'].notnull().astype(int)
loop_genes = loop_genes.reset_index(drop=True)

loop_genes.to_csv(f"{path_save}/{model}_loops_with_annotation.csv", index=False)
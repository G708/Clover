# This script queries the essentiality of human cellline data from OGEE.

# OGEE correct: 
# > highly genomically-annotated cell lines, represented by molecular features of patient tumors or alterations are gathered from multiple sources. 
# > A total of 4 datasets resulted from RNAi screen or/and CRISPR screen were collected from depmap.
# > Genome-scale CRISPR screening data are collected from sanger project score.

# As quate,
# > The essentiality of a gene is highly dependent on various factors including the genetic context, genetic background of the host and environment (15).
# > Gurumayum et al., 2021. OGEE v3: Online GEne Essentiality database with increased coverage of organisms and human cell lines. Nucleic Acids Research
# > 15. D’Elia M.A., Pereira M.P., Brown E.D. Are essential genes really essential. Trends Microbiol. 2009; 17:433–438.

# This data is contain 931 cell lines from 27 tissues.
# This data has 8 columns: genes, tissue, ess, non_ess, total, percent, dataType, source
# ess: number of essential cell lines
# non_ess: number of non-essential cell lines
# total: total number of cell lines


import pathlib
import pandas as pd

from convert import add_ids

def main(file):
	p_file = pathlib.Path(file)

	df = pd.read_csv(file, sep="\t", header=None, skiprows=1, 
							names = ["hgnc_symbol", "tissue", "ess", "non_ess", "total", "ess_percent", "dataType", "source"])
	sub_df = df[df['total']>= 10]

	mean_essentiality = sub_df.groupby('hgnc_symbol')['ess_percent'].mean()

	df_mean_essentiality = pd.DataFrame(mean_essentiality)
	result = df_mean_essentiality.reset_index()
	result.to_csv(f"{p_file.parent}/{p_file.stem}_mean.tsv", sep="\t",index=False)
	return 
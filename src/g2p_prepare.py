# Extruct human gene publication count.
# Sorse data is from NCBI
# https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz

import pathlib
import pandas as pd

from convert import add_ids

tax_id = 9606 # human

def get_human(file):
	all_g2p = pd.read_csv(file, sep="\t", header=None, skiprows=1, 
							names = ["tax_id", "GeneID", "PubMed_ID"])
	human_g2p = all_g2p[all_g2p["tax_id"]==9606]
	return human_g2p

def main(file):
	g2p = get_human(file)
	g2p_size = g2p.groupby(["GeneID"]).size().reset_index(name='N')

	p_file = pathlib.Path(file)
	g2p_size = add_ids(g2p_size,entrezgene_id = "GeneID")
	g2p_size.to_csv(f"{p_file.parent}/{p_file.stem}_human_count.txt", sep="\t",index=False)
	return

# if __name__ == '__main__':
# 	file = "/home/oba/TF_Rank_Across_Cells/git/notebooks_src/src/Tresure_hunter_20230207/tests/resources/test_g2p"
# 	main(file)
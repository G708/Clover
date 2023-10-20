import os
import argparse
import pandas as pd
import numpy as np

from data_prep import ResourceManager
from scores import glint, dowsing, treasure_hunt, ropeway
from plot_wordcloud import plot_wordcloud

import logging
logger = logging.getLogger(__name__)

"""
Parameter
"""
def input_args():
	parser = argparse.ArgumentParser(description='finding surprising DEG')
	parser.add_argument("--input", "-i",
						type=str,
						default="",
						help="input expression matrix file with FDR columns")
	parser.add_argument("--id_type",
						type=str,
						choices=['hgnc_symbol', 'ensembl_gene_id', 'entrezgene_id'],
						default=None,
						help="gene name / ID type")
	parser.add_argument("--sep",
						type=str,
						default=",",
						help="A delimiter of input")
	parser.add_argument("--gene_column",
						type=str,
						default=None,
						help="gene name / ID column name in input")	
	parser.add_argument("--fdr_column",
						type=str,
						default=None,
						help="FDR column name in input")
	parser.add_argument("--resources",
						type=str,
						default="False",
						help="Whether to download resources or not")
	parser.add_argument("--thread", "-t",
						type=int,
						default=1,
						help="thread number to rum gini_prepare.main() parallel")
	parser.add_argument("--wc_top", "-w",
						type=int,
						default=30,
						help="Rank top N gene to plot word cloud")
	parser.add_argument("--output_path","-o",
						type=str,
						default=None,
						help="An output_directory of downloaded reference files and output. Default is '~/Clover'.")	
	return parser.parse_args()

class Rank:
	"""Return Surprising DEGs scores.

	Return scores to rank surprising DEGs.

	Attributes:
		df: An expression matrix or DEG lists.
		gene_column: A column name of gene name or ID. The default is the first column of df.
		id_type: Choose from ["hgnc_symbol", "ensembl_gene_id", "entrezgene_id"].
		fdr_column: A column name of FDR. The default is the second column of df.
	"""

	def __init__(self, df, base_folder, gene_column: str = None, id_type: str = None, fdr_column: str = None):
		self.df = df
		self.id_type = id_type
		self.base_folder = base_folder

		# get resources to calculate ranking
		try:
			self.resources = self._get_resources()
		except Exception as e:
			logger.error(f"resources ERROR: {e}")

		id_choise = ['hgnc_symbol', 'ensembl_gene_id', 'entrezgene_id']
		if not (id_type in id_choise):
			raise ValueError("Invalid sim type. Expected one of: %s" % id_choise)

		else:
			self.id_type = id_type

		if gene_column is None:
			self.df.rename(columns={self.df.columns[0]: 'genename'}, inplace=True)
			self.gene_column = 'genename'
		else:
			if not (gene_column in self.df.columns):
				raise ValueError(f"column: {self.gene_column} not in input")
			else:
				self.df.rename(columns={gene_column: 'genename'}, inplace=True)
				self.gene_column = 'genename'

		if fdr_column is None:
			self.df.rename(columns={self.df.columns[1]: 'FDR'}, inplace=True)
			self.fdr_column = "FDR"

		else:
			if not (fdr_column in self.df.columns):
				raise ValueError(f"column: {fdr_column} not in input")
			else:
				self.df.rename(columns={fdr_column: 'FDR'}, inplace=True)
				self.fdr_column = "FDR"
		return

	def _get_resources(self):
		"""Loading resources and merge to df.
		"""
		resources_path = f"{self.base_folder}/resources"
		DEPrior_g2p = pd.read_csv(f"{self.base_folder}/DEPrior_gini_g2p.txt", sep="\t")		
		return DEPrior_g2p
		
	def run(self):
		"""Merge resource score and calculate ranking scores.
		
		Returns:
			df_merge: pandas.DataFrame which adds ranking columns to input df.
		"""

		df_merge = pd.merge(self.df, self.resources[["hgnc_symbol", "ensembl_gene_id", "entrezgene_id", "DE_Prior_Rank", "g2p_rank", "N", "gini_norm"]], 
			  left_on=self.gene_column, right_on=self.id_type, how="left")
		
		# Replace FDR 0 to second smallest value * 0.001
		# To avoid dividing by zero error
		unique_fdr = np.sort(df_merge["FDR"].unique())
		if 0 in unique_fdr:
			df_merge.loc[df_merge["FDR"] == 0, "FDR"] =	np.sort(df_merge["FDR"].unique())[1] * 0.001
		else:
			pass
		
		# Calculate ranking scores

		df_merge["Glint"] = glint(df_merge["DE_Prior_Rank"], df_merge["gini_norm"])
		df_merge["Dowsing"] = dowsing(df_merge["DE_Prior_Rank"], df_merge["gini_norm"], df_merge[self.fdr_column] )
		df_merge["Treasure_Hunt"] = treasure_hunt(df_merge["g2p_rank"], df_merge["Dowsing"])
		df_merge["Ropeway"] = ropeway(df_merge["g2p_rank"], df_merge["Dowsing"])

		return df_merge
	

def asc_set(score):
	asc = ['FDR', 'Glint']
	if score in asc:
		return True
	else:
		return False

def main(args):
	if args.resources == "True":
		resource = ResourceManager(thread = args.thread, base_folder = args.output_path)
		resource.download_all()
	else:
		# Find resource folder
		# if not exist, return error message
		resource = ResourceManager(thread = args.thread, base_folder = args.output_path)
		if not os.path.exists(f"{resource.base_folder}/DEPrior_gini_g2p.txt"):
			raise ValueError(f"Resource data do not exists in {resource.base_folder}: Please run with --resources True option or run `python src/data_prep.py`")

	try:
		df = pd.read_csv(args.input, sep = args.sep)
	except Exception as e:
		logger.error('input load ERROR: %s' % (e))

	rank_get = Rank(df, resource.base_folder, args.gene_column, args.id_type, args.fdr_column)
	result = rank_get.run()
	result.to_csv(f"{resource.base_folder}/rank_result.csv")
	os.makedirs(f"{resource.base_folder}/wordcloud", exist_ok=True)
	for score in ["FDR", "Glint", "Dowsing", "Treasure_Hunt", "Ropeway"]:
		sorted_result = result.sort_values(by=score, ascending=asc_set(score))
		plot_wordcloud(sorted_result[0:args.wc_top], score, f"{resource.base_folder}/wordcloud")
	return result

if __name__ == '__main__':
	args = input_args()
	main(args)
	# file = "/home/oba/TF_Rank_Across_Cells/git/notebooks_src/src/Tresure_hunter_20230207/tests/resources/test_input_dataset1.txt"
	# df = pd.read_csv(file, sep="\t")
	# df_merge = Rank(df, id_type = "hgnc_symbol").run()
	# print(df_merge)

	# file = "/home/oba/TF_Rank_Across_Cells/git/notebooks_src/src/Tresure_hunter_20230207/tests/resources/test_input_dataset2.txt"
	# df = pd.read_csv(file, sep="\t")
	# df_merge = Rank(df, id_type = "hgnc_symbol", fdr_column = "fdr").run()
	# print(df_merge)
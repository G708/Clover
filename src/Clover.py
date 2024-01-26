import os
import argparse
import pandas as pd
import numpy as np

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
						help="An output_directory of downloaded reference files and output. Default is current directory.")	
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

	def __init__(self, df, output_path, gene_column: str = None, id_type: str = None, fdr_column: str = None):
		self.df = df
		self.id_type = id_type
		self.output_path = output_path

		# output_path setting if -o option is None, work_directory = os.getcwd()
		if self.output_path is None:
			output_path = os.getcwd()
		else:	
			self.output_path = os.path.abspath(self.output_path)
			# generate output directory if exist, raise error
			if os.path.exists(self.output_path):
				raise ValueError(f"Output directory already exists: {self.output_path}")
			else:
				os.makedirs(self.output_path, exist_ok=False)

		# get resources to calculate ranking
		try:
			self.resources = self._get_resources()
		except Exception as e:
			logger.error(f"resources ERROR: {e}")

		id_choice = ['hgnc_symbol', 'ensembl_gene_id', 'entrezgene_id']
		if not (id_type in id_choice):
			raise ValueError("Invalid sim type. Expected one of: %s" % id_choice)

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
		resources_path = os.getcwd() + "/data"
		if not os.path.exists(f"{resources_path}/DEPrior_gini_g2p.txt"):
			raise ValueError(f"Resource data do not exists in {resources_path}: Please run `python src/data_prep.py`")
		else:
			DEPrior_g2p = pd.read_csv(f"{resources_path}/DEPrior_gini_g2p.txt", sep="\t")		
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

	work_directory = os.getcwd()

	try:
		df = pd.read_csv(args.input, sep = args.sep)
	except Exception as e:
		logger.error('input load ERROR: %s' % (e))

	rank_get = Rank(df, args.output_path, args.gene_column, args.id_type, args.fdr_column)
	result = rank_get.run()
	result.to_csv(f"{args.output_path}/rank_result.csv")
	os.makedirs(f"{args.output_path}/wordcloud", exist_ok=True)
	for score in ["FDR", "Glint", "Dowsing", "Treasure_Hunt", "Ropeway"]:
		sorted_result = result.sort_values(by=score, ascending=asc_set(score))
		plot_wordcloud(sorted_result[0:args.wc_top], score, f"{args.output_path}/wordcloud")
	return result

if __name__ == '__main__':
	args = input_args()
	main(args)

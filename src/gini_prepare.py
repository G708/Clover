# -*- coding: utf-8 -*-

# Get the Gini index on genes as a matrix
# Sorse RNA-seq data is from GTEx data in Human Protein Atlas
# https://www.proteinatlas.org/about/download

# output matrix has 19764 gene

import pathlib
import functools
import pandas as pd
import numpy as np
from multiprocessing import Pool

import logging
logger = logging.getLogger(__name__)

from convert import add_ids

def gini_fast(array):
	""" Calculate the Gini index

	Args:
		array: A numpy array of expression values.

	Returns:
		gini: Gini index of one gene
	"""
	# https://github.com/oliviaguest/gini
	# based on bottom eq:
	# http://www.statsdirect.com/help/generatedimages/equations/equation154.svg
	# from:
	# http://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm

	# All values are treated equally, arrays must be 1d:
	# Values cannot be negative:
	array = array.flatten()
	array -= np.amin(array)
	# Values cannot be 0:
	array += 0.0000001
	# Values must be sorted:
	array = np.sort(array)
	# Number of array elements:
	n = array.shape[0]
	# Index per array element:
	index = np.arange(1,n+1)

	gini = ((np.sum((2 * index - (n + 1)) * array)) / (n * np.sum(array)))

	return gini

class GetGini(object):
	""" Calculate Gini index

	Args: 
		file: GTEx RNA-seq file name. A column name in the matrix is as below:
				column0: Gene
				column1: Gene name
				column2: Tissue
				column3: TPM
				column4: pTPM
				column5: nTPM
	"""
	def __init__(self, file: str):

		self.file = file
		self.expression_matrix = self.read_file()
		self.type = "Tissue" # Column name to group
		self.expression_column = "nTPM" # Column name of expression value
		self.name_column = "Gene name"

	def read_file(self):
		""" Read GTEx RNA-seq expression from HPA
		
		Returns:
			pandas.DataFrame of read file.
		
		Raises:
			FileNotFoundError: An error occurs when the input file does not exist.

		"""
		try:
			df = pd.read_csv(self.file, sep="\t", 
						dtype={0: str, 1: str, 2: str, 3: float, 4: float, 5: float},
						encoding="utf-8")
			logger.info(f"Loaded: {self.file}")
			return df

		except FileNotFoundError as e:
			logger.info(f"FileNotFoundError: {e}")

	def genes(self):
		""" Extruct unique gene in matrix
		A gene must be quantitated over one tissue in the dataset.

		Returns:
			uniq_gene: list of the unique genes in the matrix which has more than one tissue data.
		"""
		logger.info(f"start gini")

		# remove genes which only have one tissue data in a matrix.
		gene_group_size = self.expression_matrix.groupby(['Gene']).size()
		gene_group = gene_group_size[gene_group_size > 1]

		# gene lists which will be returning Gini index
		uniq_gene = gene_group.index.get_level_values('Gene').unique()
		logger.info(f"Number of genes: {len(uniq_gene)}")
		return uniq_gene

	def calculate_gini(self, gene: str):
		""" Calculate gini by gene

		Args: 
			gene -> str : Name of the gene.

		Returns:
			gini: Gini index. Max is 1 - 1/N. (N is tissue size)
			g_norm: Normalized Gini index. Normalize by (N/(N-1)). Max is 1.
			expression_mean: Expression mean amoung tissue.

		"""
		subset = self.expression_matrix[self.expression_matrix['Gene'] == gene]
		expression_values = subset[self.expression_column]
		gene_name = subset[self.name_column]
		n = len(expression_values)

		gini = gini_fast(np.array(expression_values))
		g_norm = gini * (n/(n-1))
		expression_mean = expression_values.mean()

		return gene_name.unique()[0], gini, g_norm, expression_mean

def main(input_file, thread):
	get_Gini = GetGini(input_file)
	uniq_gene = get_Gini.genes()

	if thread is None:
		thread = 1
	
	with Pool(thread) as p:
		result = p.map(functools.partial(get_Gini.calculate_gini),uniq_gene)
	gene_name = [i[0] for i in result]
	gini = [i[1] for i in result]
	gini_norm = [i[2] for i in result]
	expression_mean = [i[3] for i in result]

	result = pd.DataFrame({
			'Ensembl ID': uniq_gene,
			'Gene name': gene_name,
			'gini': gini,
			"gini_norm" : gini_norm,
			'expression_mean': expression_mean})

	p_file = pathlib.Path(input_file)
	result = add_ids(result,hgnc_symbol = "Gene name", ensembl_gene_id = "Ensembl ID")
	result.to_csv(f"{p_file.parent}/{p_file.stem}_gini_norm.tsv", sep="\t",index=False)

	return 

# if __name__ == '__main__':
# 	thread = 10
# 	file = "/home/oba/TF_Rank_Across_Cells/git/notebooks_src/src/Tresure_hunter_20230207/tests/resources/test_gini_file.txt"
# 	main(file, thread)
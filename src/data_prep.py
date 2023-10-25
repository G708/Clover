# This script is used to download all the resources needed for Clover.
# gene2pubmed is downloaded from NCBI.
# rna_tissue_gtex.tsv is downloaded from Human Protein Atlas.
# DE_Prior.txt is downloaded from https://raw.githubusercontent.com/maggiecrow/DEprior/master/DE_Prior.txt.
# tissueSpecific.txt is downloaded from https://v3.ogee.info/static/files/tissueSpecific.txt.gz

# Usage:
# python3 data_prep.py -t 4 -b /path/to/base/folder
#
# -t: number of threads to run gini_prepare.main() parallel.
# -b: base folder name of downloaded reference files.
#

import os
import sys
import argparse
import logging
import pathlib
import pandas as pd
import numpy as np
from multiprocessing import Pool

import urllib.request
import gzip
import shutil
from zipfile import ZipFile
from sklearn.preprocessing import QuantileTransformer

import gini_prepare
import g2p_prepare

logger = logging.getLogger(__name__)

logging.basicConfig(level=logging.INFO,
format="%(asctime)s - %(levelname)s:%(name)s - %(message)s",
filename="download.log")

base_folder_default = "data"

class ResourceManager(object):
	"""Manage the download, caching, access, and modification of resource files.

	After Run ResourceManager.download_all(), the reference folder will created.
	The components (files) in the reference folder are as below:

		base_folder
		├── DEPrior_gini_g2p.txt
		└── resource
			├── DE_Prior.txt
			├── gene2pubmed
			├── gene2pubmed.gz
			├── gene2pubmed_human_count.txt
			├── rna_tissue_gtex_gini_norm.tsv
			├── rna_tissue_gtex.tsv
			└── rna_tissue_gtex.tsv.zip

	Attributes:
		thread: A thread number to run gini_prepare.main() parallel.
		base_folder: A base_folder name of downloaded reference files.
	"""
	def __init__(self, thread, base_folder=None):

		self.base_folder = base_folder if base_folder else \
			os.path.join(os.path.expanduser('~'), base_folder_default)
		self.resource_folder = self._get_resource_folder()
		logger.info('Using %s as resource folder.' % self.resource_folder)
		self.thread = thread

	def get_g2p(self):
		fname = os.path.join(self.resource_folder, 'gene2pubmed')
		if not os.path.exists(fname):
			url = 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz'
			download_gzip(url, fname)
		else:
			logger.warning(fname + ' already exists. Skip Download.')
		return fname

	def get_gtex_HPA(self):
		fname = os.path.join(self.resource_folder, 'rna_tissue_gtex.tsv')
		if not os.path.exists(fname):
			url = 'https://www.proteinatlas.org/download/rna_tissue_gtex.tsv.zip'
			download_zip(url, fname)
		else:
			logger.warning(fname + ' already exists. Skip Download.')
		return fname
	
	def get_DE_Prior(self):
		fname = os.path.join(self.resource_folder, 'DE_Prior.txt')
		if not os.path.exists(fname):
			url = 'https://raw.githubusercontent.com/maggiecrow/DEprior/master/DE_Prior.txt'
			download_url(url, fname)
		else:
			logger.warning(fname + ' already exists. Skip Download.')
		return fname


	def download_all(self):
		g2p_file = self.get_g2p()
		gtex_HPA_file = self.get_gtex_HPA()
		de_prior_file = self.get_DE_Prior()

		gini_preprocess = os.path.join(self.resource_folder, 'rna_tissue_gtex_gini_norm.tsv')
		g2p_preprocess = os.path.join(self.resource_folder, 'gene2pubmed_human_count.txt')

		if not os.path.exists(g2p_preprocess):
			g2p_prepare.main(g2p_file)
		else:
			logger.warning(g2p_preprocess + ' already exists. Skip preprocess.')

		if not os.path.exists(gini_preprocess):
			gini_prepare.main(gtex_HPA_file, self.thread)
		else:
			logger.warning(gini_preprocess + ' already exists. Skip preprocess.')

		if not os.path.exists(f"{self.base_folder}/DEPrior_gini_g2p.txt"):
			self.marge_all()
		else:
			logger.warning(f"{self.base_folder}/DEPrior_gini_g2p.txt" + ' already exists. Skip preprocess.')
		return 

	def marge_all(self):
		"""Loading resources and merge to df.
		"""
		g2p = pd.read_csv(f"{self.resource_folder}/gene2pubmed_human_count.txt", sep="\t", low_memory=False)
		gini = pd.read_csv(f"{self.resource_folder}/rna_tissue_gtex_gini_norm.tsv", sep="\t")
		DEPrior = pd.read_csv(f"{self.resource_folder}/DE_Prior.txt", sep="\t")

		DEPrior_g2p = pd.merge(g2p[["entrezgene_id","hgnc_symbol", "ensembl_gene_id","N"]], DEPrior[["Gene_Name", "Gene_EntrezID", "DE_Prior_Rank"]], 
			  left_on="entrezgene_id", right_on="Gene_EntrezID", how="outer")

		DEPrior_g2p = pd.merge(DEPrior_g2p, gini, 
			  on=["entrezgene_id", "ensembl_gene_id", "hgnc_symbol"], how="outer")


		DEPrior_g2p.dropna(inplace=True)

		# QuantileTransformer on gene2pubmed (uniform distribution: 0 to 1)
		# Only gene in the gene in reference matrix
		qt = QuantileTransformer()
		X = np.array(DEPrior_g2p["N"].values).reshape(-1,1)
		DEPrior_g2p["g2p_rank"] =qt.fit_transform(X)

		# QuantileTransformer on normalized Gini index (uniform distribution: 0 to 1)
		X = np.array(DEPrior_g2p["gini_norm"].values).reshape(-1,1)
		DEPrior_g2p["gini_norm_rank"] =qt.fit_transform(X)

		# Set data type
		# DEPrior_g2p = DEPrior_g2p.astype({'entrezgene_id': str})
		DEPrior_g2p.to_csv(
			f"{self.base_folder}/DEPrior_gini_g2p.txt", sep="\t",index=False,
			columns = [
				'hgnc_symbol', 'entrezgene_id', 'ensembl_gene_id',
				'DE_Prior_Rank', 'N', 'g2p_rank',
				'gini', 'gini_norm', 'gini_norm_rank'
			]
		)
		
		return DEPrior_g2p

	def _get_resource_folder(self):
		resource_dir = os.path.join(self.base_folder, 'resources')
		if not os.path.isdir(resource_dir):
			try:
				os.makedirs(resource_dir)
			except Exception:
				logger.warning(resource_dir + ' already exists.')
		return resource_dir
 

def download_url(url, fname):
	try:
		logger.info('Downloading %s into %s' % (url, fname))
		urllib.request.urlretrieve(url, fname)
	except Exception as e:
		logger.error('Download ERROR: %s' % (e))


def download_gzip(url, fname):
	try:
		logger.info('Downloading %s and extracting into %s' % (url, fname))
		gz_file = fname + ".gz"
		urllib.request.urlretrieve(url, gz_file)
		with gzip.open(gz_file, 'rb') as fin:
			with open(fname, 'wb') as fout:
				shutil.copyfileobj(fin, fout)
	except Exception as e:
		logger.error('Download ERROR: %s' % (e))
		
def download_zip(url, fname):
	try:
		logger.info('Downloading %s and extracting into %s' % (url, fname))
		zfile = fname + ".zip"
		urllib.request.urlretrieve(url, zfile)
		with ZipFile(zfile, "r") as f:
			f.extractall(path=os.path.dirname(fname))
	except Exception as e:
		logger.error('Download ERROR: %s' % (e))

if __name__ == '__main__':
	# Download all the resources if this script is run directly
	parser = argparse.ArgumentParser(description='Download resources')
	parser.add_argument("--thread", "-t",
                    type=int,
                    default=1,
                    help="thread number to rum gini_prepare.main() parallel")
	parser.add_argument("--output_path","-o",
                    type=str,
                    default=None,
                    help="A output_directory of downloaded reference files and output. Default is '~/Clover'.")
	args = parser.parse_args()
	
	ResourceManager(args.thread, args.output_path).download_all()

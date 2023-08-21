# This script queries the essentiality of human cellline data from OGEE.

import pathlib
import pandas as pd

from convert import add_ids



def main(file):
	df = pd.read_csv(file, sep="\t", header=None, skiprows=1, 
							names = ["hgnc_symbol", "tissue", "ess", "non_ess", "total", "percent", "dataType", "source"])
	sub_df = df[df['total']>= 10]


	
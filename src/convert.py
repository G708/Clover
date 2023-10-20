from biomart import BiomartServer
import pandas as pd

import logging
logger = logging.getLogger(__name__)

def get_ensembl_mappings():
	"""Connect ensembl biomart through API.
	
	Returns:
		to_genesymbol: Dictionary which keys are ensembl_gene_id or entrezgene_id and value is hgnc_symbol.
		to_entrezgene: Dictionary which keys are ensembl_gene_id or hgnc_symbol and value is entrezgene_id.
		to_ensembl: Dictionary which keys are hgnc_symbol or entrezgene_id and value is ensembl_gene_id.
	"""
	# Set up a connection to the server. Change the mirror site to your location.
	server = BiomartServer( "http://asia.ensembl.org/biomart" )
	# server = BiomartServer( "http://www.ensembl.org/biomart" )
	mart = server.datasets['hsapiens_gene_ensembl']
	
	attributes = ['hgnc_symbol', 'ensembl_gene_id', "entrezgene_id"]
	response = mart.search({'attributes': attributes}) 
	data = response.raw.data.decode('ascii')

	to_genesymbol = {}
	to_entrezgene = {}
	to_ensembl = {}
	# Store the data in a dict
	for line in data.splitlines():
		line = line.split('\t')
		# The entries are in the same order as in the `attributes` variable
		gene_symbol = line[0]
		ensembl_gene = line[1]
		entrezgene_id = line[2]

		# Some of these keys may be an empty string.
		to_genesymbol[ensembl_gene] = gene_symbol
		to_genesymbol[entrezgene_id] = gene_symbol

		to_entrezgene[ensembl_gene] = entrezgene_id
		to_entrezgene[gene_symbol] = entrezgene_id

		to_ensembl[entrezgene_id] = ensembl_gene
		to_ensembl[gene_symbol] = ensembl_gene

	return to_genesymbol, to_entrezgene, to_ensembl

def add_ids(df, hgnc_symbol = "hgnc_symbol", ensembl_gene_id = "ensembl_gene_id", entrezgene_id = "entrezgene_id"):
	"""Add hgnc_symbol, ensembl_gene_id, or entrezgene_id.

	This will add hgnc_symbol, ensembl_gene_id, or entrezgene_id to df if there is no corresponding column.
	Also, unify column names.

	Args:
		df: pandas.DataFrame.
		hgnc_symbol: corresponding column name to hgnc_symbol (gene name).
		ensembl_gene_id: corresponding column name to ensembl_gene_id, start from ENSG**********.
		entrezgene_id: corresponding column name to entrezgene_id (only number).
	
	Returns:
		df_modify: pandas.DataFrame which added an id and symbol, also, changes column names.
	"""
	attributes = {
		hgnc_symbol: "hgnc_symbol", 
		ensembl_gene_id: "ensembl_gene_id", 
		entrezgene_id: "entrezgene_id"}
    
	df_modify = df.rename(columns=attributes)
	df_modify

	to_genesymbol, to_entrezgene, to_ensembl = get_ensembl_mappings()

	ids_in = set(df_modify.columns) & set(attributes.values())
	ids_not_in = set(attributes.values()) - ids_in

	for attr in ids_not_in:
		ref = list(ids_in)[0]
        
		if attr == hgnc_symbol: # convert from ensembl_gene_id or entrezgene_id to hgnc_symbol
			logger.info(f"convert hgnc_symbol from {ref}")
			df_modify[attr] = df_modify[ref].astype(str).map(to_genesymbol)
            
		elif attr == ensembl_gene_id: # convert from hgnc_symbol or entrezgene_id to ensembl_gene_id
			logger.info(f"convert ensembl_gene_id from {ref}")
			df_modify[attr] = df_modify[ref].astype(str).map(to_ensembl)
            
		elif attr == entrezgene_id: # convert from hgnc_symbol or ensembl_gene_id to entrezgene_id
			logger.info(f"convert entrezgene_id from {ref}")
			df_modify[attr] = df_modify[ref].astype(str).map(to_entrezgene)
            
		else:
			logger.error("convert error")
	return df_modify

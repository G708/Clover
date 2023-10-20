# Clover

## Command line usage

```bash
usage: Clover.py [-h] [--input INPUT] [--id_type {hgnc_symbol,ensembl_gene_id,entrezgene_id}]
				[--sep SEP] [--gene_column GENE_COLUMN]
				[--fdr_column FDR_COLUMN] [--resources RESOURCES]
				[--thread THREAD] [--wc_top WC_TOP] [--output_path OUTPUT_PATH]

```

finding surprising DEG

## optional arguments:

### `-h`, `--help`

show this help message and exit.

### `--input` / `-i` INPUT

[**Required Arguments**] Input file of the DEG list with FDR columns.

### `--id_type`

[**Required Arguments**] gene name / ID type to convert ID and merge files. Choose from;{hgnc_symbol, ensembl_gene_id, entrezgene_id}.

### `--sep` SEP

A delimiter of input. Default is `,`.

### `--gene_column` GENE_COLUMN

Gene name/ID column name in the input file. If this is not specified, Clover will recognize the first column as the gene name/ID.

### `--fdr_column` FDR_COLUMN

FDR column name in input. If this is not specified, Clover will recognize the second column as the FDR.

### `--resources` RESOURCES

Whether to download resources or not. Default is `False`.

### `--thread` / `-t` THREAD

Thread number to rum `gini_prepare.main()` parallel. Default is 1.

### `--wc_top` / `-w` WC_TOP

Rank top N gene to plot word cloud. Default is 30.

### `--output_path` / `-o` OUTPUT_PATH

A output_directory of downloaded reference files and output. Default is `~/Clover`.

# Clover_resources

## ResourceManager()

Manage the download, caching, access, modification of resource files.

After runing `ResourceManager.download_all()`, reference folder will created.
Components (files) in reference folder is as below:

```bash
output_path
├── DEPrior_gini_g2p.txt
└── resources
	├── DE_Prior.txt
	├── gene2pubmed
	├── gene2pubmed.gz
	├── gene2pubmed_human_count.txt
	├── rna_tissue_gtex_gini_norm.tsv
	├── rna_tissue_gtex.tsv
	└── rna_tissue_gtex.tsv.zip
```

### Attributes:

- thread: A thread number to run `gini_prepare.main()` parallel.
- base_folder: A base_folder name of downloaded reference files

# plot_wordcloud

This function will plot the wordcloud plot based on the ranking and the score.

```python
plot_wordcloud(df, score, output_path, figsize=(10, 3), colormap=plt.get_cmap("viridis"))
```

### Attributes:

- df: score matrix culculated by Clover
- score: Score name which is column of the df
- output_path: output path name
- figsize: Figure size. Default is (10, 3).
- colormap: colormap of the font color. Default is plt.get_cmap("viridis").

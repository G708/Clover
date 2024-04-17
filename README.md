# Clover

## 0. ChangeLog
Please see [ChangeLog.md](ChangeLog.md)

## 1. Overview
`Clover` is a tool for prioritising rare genes in a differentially expressed gene (DEG) list.

## 2. Installation

### 2.1. Requirements
* Python 3.8 or higher
Other requirements are listed in [requirements.txt](requirements.txt).


To run `Clover`, clone this repository and move to the repository in your terminal / Command Prompt.
```bash
git clone https://github.com/G708/Clover.git
cd Clover
```

The project relies on several Python packages. You can install these dependencies like this:
```bash
pip install -r requirements.txt
```

Download the data with:
```bash
python src/data_prep.py
```

## 3. Quick usage
This tool can use as command line tool and interactive dashboard application.


### 3.1 Command line tool
To run Clover as command line tool, run this command:

```bash
python src/Clover.py -i test_imput_DEG.csv --id_type hgnc_symbol -o test_output
```

#### Output
The output files are as below:
  
  ```bash
  Clover
  └──test_output
    ├──  rank_result.csv
    └── wordcloud
      ├── FDR.png
      ├── Glint.png
      ├── Dowsing.png
      ├── Treasure_Hunt.png
      └── Ropeway.png
  ```

#### 3.1.1 Essential options
To run command line tool, you need to specify input file with `-i` option.
This input file should be a csv file with FDR columns. For full options, please refer to [ function page](docs/functions.md).

##### `-i`, `--input`
Input expression matrix file with FDR columns. Default delimiter of the file is `,` , but user can specify by the additional option `--sep`.

##### `--id_type`
gene name / ID type in input file to convert ID and merge the columns. Choose from `{hgnc_symbol, ensembl_gene_id, entrezgene_id}`. Default is `hgnc_symbol`.

### 3.2. Interactive Dash tool

For interactive dashboard application, we implemented clover using [Dash](https://dash.plotly.com/).

Please check the [https://github.com/G708/clover-dash](https://github.com/G708/clover-dash) repository for more details.


## 4. Feature Data sources

Download the data with:
```bash
python src/data_prep.py
```

### Data sources information

* RNA-seq data from Genotype-Tissue Expression GTEx project v8 (THE GTEX CONSORTIUM, 2020) data in Human protein Atlas ([https://www.proteinatlas.org/download/rna_tissue_gtex.tsv.zip](https://www.proteinatlas.org/download/rna_tissue_gtex.tsv.zip)).
* Publication data for each gene is publicly available at [https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz](https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz).
* Score of DE prior is downloaded from
  > Crow,M. et al. (2019) Predictability of human differential gene expression. Proceedings of the National Academy of Sciences, 116, 6491–6500. [https://doi.org/10.1073/pnas.1802973116](https://doi.org/10.1073/pnas.1802973116)

## 5. Citation

> Oba GM, Nakato R. Clover: An unbiased method for prioritizing differentially expressed genes using a data-driven approach. Genes Cells. 2024 Apr 11. doi: [10.1111/gtc.13119](https://doi.org/10.1111/gtc.13119). Epub ahead of print. PMID: 38602264.

## 6. License
This project is licensed under the terms of the MIT license.

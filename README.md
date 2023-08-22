# Clover-vis
`Clover` tool to mining differentially expressed genes (DEGs) in intaractive visualization tool.
This tool helps to summarise gene feature of your intrust.

## Installation

To run this app locally, clone this repository and open this app folder in your terminal/Command Prompt.
```bash
git clone https://github.com/G708/Clover.git
cd Clover
```

The project relies on several Python packages including Dash, Dash Bootstrap Components, Plotly, and Pandas.
You can install these dependencies like this:
```bash
pip install -r requirements.txt
```

Download the data with:
```bash
python src/data_prep.py
```


Run this app locally with:
```bash
python app.py
```

Open a browser at http://127.0.0.1:8050

## Quick usase




## Feature Data sources

* RNA-seq data from Genotype-Tissue Expression GTEx project v8 (THE GTEX CONSORTIUM, 2020) data in Human protein Atlas ([https://www.proteinatlas.org/download/rna_tissue_gtex.tsv.zip](https://www.proteinatlas.org/download/rna_tissue_gtex.tsv.zip)).
* Publication data for each gene is publicly available at [https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz](https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz).
* Score of DE prior is downloaded from
  > Crow,M. et al. (2019) Predictability of human differential gene expression. Proceedings of the National Academy of Sciences, 116, 6491–6500. [https://doi.org/10.1073/pnas.1802973116](https://doi.org/10.1073/pnas.1802973116)
  >
* Human cell line essentiality data is downloaded from OGEE (https://v3.ogee.info/#/home)
  > S. Gurumayum *, P. Jiang *, X. Hao, T.L Campos, N.D. Young, P.K. Korhonen, R.B. Gasser, P. Bork, X.M. Zhao, L.J. He, Chen, W.-H. (2020). "OGEE v3: Online GEne Essentiality database with increased coverage of organisms and human cell lines. Nucleic Acids Research, gkaa884, [doi: 10.1093/nar/gkaa884 ](https://doi.org/10.1093/nar/gkaa884)
  >

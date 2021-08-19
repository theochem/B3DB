# About *B3DB*

In this repo, we present a large benchmark dataset, Blood-Brain Barrier Database (B3DB), complied
from 50 published resources (as summaried at
[raw_data/raw_data_summary.tsv](raw_data/raw_data_summary.tsv)) and categorized based on
experimental uncertainty. A subset of the
molecules in B3DB has numerical `logBB` values (1058 compounds), while the whole dataset
has categorical (BBB+ or BBB-) BBB permeability labels (7807). Some physicochemical properties
of the molecules are also provided.

## Citation

Please use the following citation in any publication using Procrustes library:

> **"B3DB: A Curated Database of Blood-Brain Barrier Permeability and Chemical Descriptors for a
Diverse Set of Compounds", F. Meng, et al.**

*To be updated once the publication is out.*

## Features of *B3DB*

1. The largest dataset with numerical and categorical values for Blood-Brain Barrier small molecules
    (to the best of our knowledge as of 2021/Feb/25).

2. Inclusion of sterochemistry information with isomeric SMILES with chiral specifications if
    available. Otherwise, canonical SMILES are used.

3. Characterization of uncertainty of experimental measurements by grouping the collected molecular
    data records.

4. Extended datasets for numerical and categorical data with precomputed physicochemical properties
    using [mordred](https://github.com/mordred-descriptor/mordred).

## Usage

There are two types of dataset in [B3DB](B3DB), [regression data](B3DB/B3DB_regression.tsv)
and [classification data](B3DB/B3DB_classification.tsv) and they can just simply load with *pandas*
library. For example

```python
import pandas as pd

# load regression data
regression_data = pd.read_csv("B3DB/B3DB_regression.tsv",
                              sep="\t")

# load classification data
classification_data = pd.read_csv("B3DB/B3DB_classification.tsv",
                                  sep="\t")
```

## Working environment setting up

All the calculations were performed in a Python 3.7.9 virtual environment created with Conda in
CentOS Linux release 7.9.2009 includes Python packages,

- ChEMBL_Structure_Pipeline==1.0.0, https://github.com/chembl/ChEMBL_Structure_Pipeline/
- RDKit==2020.09.1, https://www.rdkit.org/
- openeye-toolkit===2020.2.0, https://docs.eyesopen.com/toolkits/python/index.html/
- mordred==1.1.2, https://github.com/mordred-descriptor/mordred/ (required networkx==2.3.0)
- numpy==1.19.2, https://numpy.org/
- pandas==1.2.1, https://pandas.pydata.org/
- pubchempy==1.0.4, https://github.com/mcs07/PubChemPy/
- PyTDC==0.1.5, https://github.com/mims-harvard/TDC/
- SciPy==1.5.2, https://www.scipy.org/
- tabula-py==2.2.0, https://pypi.org/project/tabula-py/

We will create a virtual environment named *bbb_data* with `Python 3.7.9` first,
```bash
conda create bbb_py37 python=3.7.9
```
Given that `RDKit`, `ChEMBL_Structure_Pipeline` are not available in PyPI and we will install
them with `conda`,

```bash
# activate virtual environment
conda activate bbb_py37

conda install -c rdkit rdkit=2020.09.1.0
conda install -c conda-forge chembl_structure_pipeline=1.0.0
# https://docs.eyesopen.com/toolkits/python/quickstart-python/linuxosx.html
conda install -c openeye openeye-toolkits=2020.2.0
```
Then we can install the requirements in [requirements.txt](requirements.txt) with
```bash
pip install -r requirements.txt
```

An easier way is to run the follow script with `bash`,

```bash
#!/bin/bash

# create virtual environment
conda create bbb_py37 python=3.7.9
# activate virtual environment
conda activate bbb_py37

# install required packages
conda install -c rdkit rdkit=2020.09.1.0
conda install -c conda-forge chembl_structure_pipeline=1.0.0
pip install -r requirements.txt
```

`ALOGPS` version 2.1 can be accessed at http://www.vcclab.org/lab/alogps/.

The materials and data under this repo are under
[CC-BY-4.0 Licence](https://creativecommons.org/licenses/by/4.0/legalcode).

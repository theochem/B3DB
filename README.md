# About *B3DB*

In this repo, we present a large benchmark dataset, [Blood-Brain Barrier Database (B3DB)](https://www.nature.com/articles/s41597-021-01069-5), compiled
from 50 published resources (as summarized at
[raw_data/raw_data_summary.tsv](raw_data/raw_data_summary.tsv)) and categorized based on
the consistency between different references/measurements. This dataset was [published in Scientific Data](https://www.nature.com/articles/s41597-021-01069-5) and this repository is occassionally uploaded with new data. Scientists who would like to contribute data should contact the database's maintainers (e.g., by creating a new Issue in this database).

A subset of the
molecules in B3DB has numerical `logBB` values (1058 compounds), while the whole dataset
has categorical (BBB+ or BBB-) BBB permeability labels (7807 compounds). Some physicochemical properties
of the molecules are also provided. 

## Citation

Please use the following citation in any publication using our *B3DB* dataset:

```md
@article{Meng_A_curated_diverse_2021,
author = {Meng, Fanwang and Xi, Yang and Huang, Jinfeng and Ayers, Paul W.},
doi = {10.1038/s41597-021-01069-5},
journal = {Scientific Data},
number = {289},
title = {A curated diverse molecular database of blood-brain barrier permeability with chemical descriptors},
volume = {8},
year = {2021},
url = {https://www.nature.com/articles/s41597-021-01069-5},
publisher = {Springer Nature}
}
```

## Features of *B3DB*

1. The largest dataset with numerical and categorical values for Blood-Brain Barrier small molecules
    (to the best of our knowledge, as of February 25, 2021).

2. Inclusion of stereochemistry information with isomeric SMILES with chiral specifications if
    available. Otherwise, canonical SMILES are used.

3. Characterization of uncertainty of experimental measurements by grouping the collected molecular
    data records.

4. Extended datasets for numerical and categorical data with precomputed physicochemical properties
    using [mordred](https://github.com/mordred-descriptor/mordred).

## Usage

There are two types of dataset in [B3DB](B3DB), [regression data](B3DB/B3DB_regression.tsv)
and [classification data](B3DB/B3DB_classification.tsv) and they can be loaded simply using *pandas*. For example

```python
import pandas as pd

# load regression dataset
regression_data = pd.read_csv("B3DB/B3DB_regression.tsv",
                              sep="\t")

# load classification dataset
classification_data = pd.read_csv("B3DB/B3DB_classification.tsv",
                                  sep="\t")

# load extended regression dataset
regression_data_extended = pd.read_csv("B3DB/B3DB_regression_extended.tsv.gz", 
                                       sep="\t", compression="gzip")

# load extended classification dataset
classification_data_extended = pd.read_csv("B3DB/B3DB_classification_extended.tsv.gz",
                                           sep="\t", compression="gzip")

```

We also have three examples to show how to use our dataset, 
[numerical_data_analysis.ipynb](notebooks/numerical_data_analysis.ipynb), 
[PCA_projection_fingerprint.ipynb](notebooks/PCA_projection_fingerprint.ipynb) and 
[PCA_projection_descriptors.ipynb](notebooks/PCA_projection_descriptors.ipynb). 
[PCA_projection_descriptors.ipynb](notebooks/PCA_projection_descriptors.ipynb) uses precomputed 
chemical descriptors for visualization of chemical space of `B3DB`, and can be used directly
using *MyBinder*,
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/theochem/B3DB/main?filepath=notebooks%2FPCA_projection_descriptors.ipynb).
Due to the difficulty of installing `RDKit` in *MyBinder*, only `PCA_projection_descriptors.
ipynb` is set up in *MyBinder*.

## Working environment setting up

All the calculations were performed in a Python 3.7.9 virtual environment created with Conda in
CentOS Linux release 7.9.2009. The Conda environment includes the following Python packages,

- ChEMBL_Structure_Pipeline==1.0.0, https://github.com/chembl/ChEMBL_Structure_Pipeline/
- RDKit==2020.09.1, https://www.rdkit.org/
- openeye-toolkit==2020.2.0, https://docs.eyesopen.com/toolkits/python/index.html/
- mordred==1.1.2, https://github.com/mordred-descriptor/mordred/ (required networkx==2.3.0)
- numpy==1.19.2, https://numpy.org/
- pandas==1.2.1, https://pandas.pydata.org/
- pubchempy==1.0.4, https://github.com/mcs07/PubChemPy/
- PyTDC==0.1.5, https://github.com/mims-harvard/TDC/
- SciPy==1.5.2, https://www.scipy.org/
- tabula-py==2.2.0, https://pypi.org/project/tabula-py/

To creat a virtual environment named *bbb_data* with `Python 3.7.9` to this specification, first,
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
# https://docs.eyesopen.com/toolkits/python/quickstart-python/linuxosx.html
conda install -c openeye openeye-toolkits=2020.2.0

pip install -r requirements.txt
```

`ALOGPS` version 2.1 can be accessed at http://www.vcclab.org/lab/alogps/.

The materials and data under this repo are distributed under the
[CC0 Licence](http://creativecommons.org/publicdomain/zero/1.0/).

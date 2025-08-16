# About *B3DB*

In this repo, we present a large benchmark dataset, [Blood-Brain Barrier Database (B3DB)](https://www.nature.com/articles/s41597-021-01069-5), compiled
from 50 published resources (as summarized at
[raw_data/raw_data_summary.tsv](raw_data/raw_data_summary.tsv)) and categorized based on
the consistency between different experimental references/measurements. This dataset was [published in Scientific Data](https://www.nature.com/articles/s41597-021-01069-5) and this repository is occasionally uploaded with new experimental data. Scientists who would like to contribute data should contact the database's maintainers (e.g., by creating a new Issue in this database).

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

### Via PyPI

The B3DB dataset is avaliable at [PyPI](https://pypi.org/project/qc-B3DB/). One can install it using pip:

```bash
pip install qc-B3DB
```

Then load the data (dictionary of `pandas` dataframe) with the following code snippet:

```python

from B3DB import B3DB_DATA_DICT

# access the data via dictionary keys
# 'B3DB_regression'
# 'B3DB_regression_extended'
# 'B3DB_classification'
# 'B3DB_classification_extended'

```

### Manually Download the Data

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

### Examples in Jupyter Notebooks

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

## Data Curation

Detailed procedures for data curation can be found in [data curation section](data_curation/) in this repository.

The materials and data under this repo are distributed under the
[CC0 Licence](http://creativecommons.org/publicdomain/zero/1.0/).

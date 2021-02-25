# About *B3DB*

In this repo, we present a large benchmark dataset, Blood-Brain Barrier Database (B3DB), complied
from 50 published resources and categorized based on experimental uncertainty. A subset of the
molecules in B3DB has numerical `logBB` values (1058 compounds), while the whole dataset
has categorical (BBB+ or BBB-) BBB permeability labels (7807). Some physicochemical properties
of the molecules are also provided, .

## Citation

Please use the following citation in any publication using Procrustes library:

> **"B3DB: A Curated Database of Blood-Brain Barrier Permeability and Chemical Descriptors for a
Diverse Set of Compounds", F. Meng, et al.**

*To be updated once the publications is out.*

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
and [classification data](B3DB/B3DB_classification.tsv) and they can just simply load with *pands*
library. For example

```python
import pandas as pd

regression_data = pd.read_csv("B3DB/B3DB_regression.tsv",
                              sep="\t")

```

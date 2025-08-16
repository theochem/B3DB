# Data Curation Process of B3DB

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
- SciPy==1.10.0, https://www.scipy.org/
- tabula-py==2.2.0, https://pypi.org/project/tabula-py/

To creat a virtual environment named *bbb_data* with `Python 3.7.9` to this specification, first,
```bash
conda create bbb_py37 python=3.7.9
```
Given that `RDKit`, `ChEMBL_Structure_Pipeline` are not available in PyPI and we will install
them with `conda`,

```bash
# activate a virtual environment
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

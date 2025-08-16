import argparse
import os

import pandas as pd
# from rdkit import Chem
import pubchempy as pcp
from urllib.request import urlopen

__version__ = "0.0.4"


def add_CID_again(input_excel="2_bbb_all_complete_no_spaces.xlsx",
                  output_excel="2_bbb_all_complete_CID.xlsx"):
    """This function adds CID information to excel file again."""
    df = pd.read_excel(input_excel)

    counter = 0
    for idx, row in df.iterrows():
        if pd.isna(row["CID"]):
            try:
                compound = pcp.get_compounds(identifier=row["smiles"],  namespace="smiles")
                if len(compound) >= 2:
                    df.loc[idx, "comments"] = row["comment"] + \
                                              "More than one CID available  for this compound."
                    # add cid information
                df.loc[idx, "CID"] = compound[0].cid
                if compound[0].cid is not None:
                    counter += 1

            except:
                pass
    print(counter)

    df.to_excel(output_excel, index=None, engine="openpyxl")


def add_isomeric(input_excel,
                 out_excel=None):
    """Add isomeric smiles."""
    # df = pd.read_excel("2_bbb_all_complete.xls")
    df = pd.read_excel(input_excel)

    if out_excel is None:
        out_excel = input_excel.strip(".xlsx") + "_out.xlsx"

    if "isomeric_smiles" not in df.columns.to_list():
        df["isomeric_smiles"] = ""
    if "iupac_name" not in df.columns.tolist():
        df["iupac_name"] = ""
    if "smiles_result" not in df.columns.tolist():
        df["smiles_result"] = ""

    for idx, row in df.iterrows():
        try:
            if not pd.isna(row["CID"]):
                # compound = pcp.Compound.from_cid(int(row["CID"]))
                url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/property/isomericsmiles/txt".format(int(row["CID"]))
                isomeric_smiles = urlopen(url).read().decode('utf8').strip()
                # add isomeric smiles
                df.loc[idx, "isomeric_smiles"] = isomeric_smiles
                # add iupac names
                url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/property/IUPACName/txt".format(int(row["CID"]))
                df.loc[idx, "iupac_name"] = urlopen(url).read().decode('utf8').strip()

                # if CID is none, isomeric smiles is none
                df.loc[idx, "smiles_result"] = isomeric_smiles
            else:
                df.loc[idx, "smiles_result"] = row["smiles"]
        except:
            pass


    if out_excel is not None:
        df.to_excel(out_excel, index=None, engine="openpyxl")
    return df


def add_smiles_result(input_excel="2_bbb_all_complete_CID_out.xlsx",
                      out_excel="2_bbb_all_complete_CID_out_smiles.xlsx"):
    df = pd.read_excel(input_excel)
    for idx, row in df.iterrows():
        if pd.isna(row["smiles_result"]):
            df.loc[idx, "smiles_result"] = row["smiles"]

    df.to_excel(out_excel, index=None, engine="openpyxl")


if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
    # # Add the positional parameter
    # parser.add_argument('--input_excel', type=str)

    # # Parse the arguments
    # args = parser.parse_args()

    add_CID_again(input_excel="2_bbb_all_complete_no_spaces.xlsx",
                  output_excel="2_bbb_all_complete_CID.xlsx")
    add_isomeric(input_excel="2_bbb_all_complete_CID.xlsx")

    add_smiles_result(input_excel="2_bbb_all_complete_CID_out.xlsx",
                      out_excel="2_bbb_all_complete_CID_out_smiles.xlsx")

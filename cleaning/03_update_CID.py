import argparse
import os

import pandas as pd
# from rdkit import Chem
import pubchempy as pcp

__version__ = "0.0.5"


def update_CID(input_excel,
               output_excel):
    """This function updates CID information."""
    df = pd.read_excel(input_excel)

    counter = 0
    for idx, row in df.iterrows():
            try:
                compound = pcp.get_compounds(identifier=row["smiles_fixed_rdkit"], namespace="smiles")
                if len(compound) >= 2:
                    df.loc[idx, "comments"] = row["comment"] + \
                                              "More than one CID available  for this compound."
                # add cid information
                df.loc[idx, "CID"] = compound[0].cid
                if compound[0].cid is not None:
                    counter += 1

            except:
                print(row["NO."])

    print(counter)

    df.to_excel(output_excel, index=None, engine="openpyxl")


if __name__ == "__main__":
    update_CID(input_excel="2_bbb_all_complete_CID_out_smiles_fixed.xlsx",
               output_excel="2_bbb_all_complete_CID_out_smiles_fixed_updated.xlsx")

# combine files together
# 1. get a set of non-duplicated items and duplicated items by checking inchi values
# where the duplicated items have more than one record

# 2. sort the excel by inchi, logBB, compound_name
# 3.

# clean data:
# 1. standardize
# 2. remove duplicates
# 3. wrap into a good shape: add isomeric smiles; use canonical smiles;

# cid can be used to remove duplicates

import os

import pandas as pd
from rdkit import Chem
import pubchempy as pcp


def combine_excels(wdir="/home/legend/Desktop/work/BBB_database/BBB_database/references/data_src",
                   out_excel="0_bbb_all.xls"):
    """Combine all the excel files."""
    excel_paths = [os.path.join(wdir, f, "data_formatted_done.xls") for f in os.listdir(wdir)
                   if f.startswith("R")]
    # sort paths by file names
    excel_paths.sort(key=lambda x: int(x.split("/")[-2][1:]))
    frames = [pd.read_excel(f) for f in excel_paths]
    result = pd.concat(frames)
    result.to_excel(out_excel, index=None)
    return result


def remove_nan(df,
               out_excel="2_bbb_all_complete.xls"):
    # remove item with missing smiles
    df_complete = df[df.smiles.notna()]
    df_complete.to_excel(out_excel)
    return


def update_inchi(df,
                 out_excel="bbb_all_with_inchi.xls"):
    """Update inchi values."""

    # compute the inchi values
    for idx, row in df.iterrows():
        try:
            mol = Chem.MolFromSmiles(row["smiles"])
            if mol is not None:
                inchi_value = Chem.MolToInchi(mol)
                df.loc[idx, "Inchi"] = inchi_value
        except:
            df.loc[idx, "comments"] = row["comments"] + "|error of generating inchi"

    df.to_excel(out_excel, index=None)
    return df


def split_regression_classification(input_excel="2_bbb_all_complete_no_space.xls",
                                    regression_file="regression_all.xlsx",
                                    classification_file="classification_all.xlsx"):
    df = pd.read_excel(input_excel)
    df_regression = df[~df["logBB"].isnull()]
    df_regression.to_excel(regression_file, engine="openpyxl")

    # df_classification = df[df["logBB"].isnull()]
    # df_classification.to_excel(classification_file, engine="openpyxl")

    # return df_regression, df_classification
    return df_regression


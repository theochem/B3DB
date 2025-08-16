# preprocessing of all the data in each folder
from pathlib import Path

import pandas as pd
import pubchempy as pcp


# when smiles available, use smiles to get CID, CAS, (Inchi)
# when smiles not available, use name to get smiles, CID, CAS, (Inchi)
# # when multiple results available for one smiles, such as 'quinidine' has multiple CID


def named_2smi_cid(input_excel="data_formatted.xls",
                   out_excel=None,
                   add_reference=True):
    """
    Compound name resolver.

    Parameters
    ----------
    input_excel : str, optional
        Input excel file name., by default "data_formatted.xls"
    out_excel : [type], optional
        Output excel file name, by default None
    """
    if out_excel is None:
        out_excel = "data_formatted_smiles.xls"

    df = pd.read_excel(input_excel)
    # add a new column for comments
    if "comments" not in df.columns.tolist():
        df["comments"] = ""
    # add reference information if required
    if add_reference and "reference" not in df.columns.tolist():
        path = Path.cwd()
        ref_id = path.name.split(" ")[0]
        df["reference"] = ref_id

    for idx, row in df.iterrows():
        compound_name = row["compound_name"].lower()
        if compound_name is not None:
            # smiles
            try:
                # smi = cirpy.resolve(compound_name, "smiles")
                compound = pcp.get_compounds(compound_name, "name", record_type="2d")

                # if return more than one CID, will double check if the first one is good,
                # else, write a comment
                if len(compound) >= 2:
                    pubchem_names = compound[0].synonyms
                    pubchem_names = [name.lower() for name in pubchem_names]
                    if compound_name not in pubchem_names:
                        df.loc[idx, "comments"] = "More than one CID available for {}|{}".format(
                            compound_name, compound)
                smi = compound[0].canonical_smiles
                # even if smi is None, the following line still works
                df.loc[idx, "smiles"] = smi
                # CID
                df.loc[idx, "CID"] = compound[0].cid
            except:
                print("can not resolve smiles for {}: {}".format(idx + 2, compound_name))
            # cas
            # cas_number = cirpy.resolve(compound_name, "smiles")
        else:
            print("compound name is None for {}".format(idx + 2))
    df.to_excel(out_excel, index=None)


# data_folders = [f for f in next(os.walk("."))[1] if not
#                 f.startswith(("Understanding", "not", "preprocessing"))]
# data_folders.sort(key= lambda x: int(x.split(" ")[0].strip("R")))

# for f in data_folders:
#     os.chdir(f)
#     add_inchi()
#     os.chdir("..")


def cid2smi(input_excel="data_formatted.xls",
            out_excel=None,
            add_reference=True):
    """Add information based on CID of pubchem."""
    if out_excel is None:
        out_excel = "data_formatted_smiles.xls"

    df = pd.read_excel(input_excel)
    # add a new column for comments
    if "comments" not in df.columns.tolist():
        df["comments"] = ""
    # add reference information if required
    if add_reference and "reference" not in df.columns.tolist():
        path = Path.cwd()
        ref_id = path.name.split(" ")[0]
        df["reference"] = ref_id

    for idx, row in df.iterrows():
        cid = row["CID"]
        if cid is not None:
            compound = pcp.Compound.from_cid(cid)
            smi = compound.canonical_smiles
        df.loc[idx, "smiles"] = smi

    df.to_excel(out_excel, index=None)


def smi2cid_name(input_excel="data_formatted_inchi.xls",
                 out_excel="data_formatted_done.xls",
                 add_reference=True,
                 add_new_names=True):
    """Get CID from smiles."""
    # if out_excel is None:
    #     out_excel = "data_formatted_smiles.xls"

    df = pd.read_excel(input_excel)

    # add a new column for comments
    if "comments" not in df.columns.tolist():
        df["comments"] = ""
    # add reference information if required
    if add_reference and "reference" not in df.columns.tolist():
        path = Path.cwd()
        ref_id = path.name.split(" ")[0]
        df["reference"] = ref_id

    if add_new_names and "new_name" not in df.columns.tolist():
        df["new_name"] = ""

    for idx, row in df.iterrows():
        smi = row["smiles"]
        if pd.isna(smi):
            print("Smiles is empty for {}".format(idx + 2))
        else:
            # if not pd.isnull(row["CID"]) and row["CID"] != "":
            if not pd.isnull(row["CID"]):
                try:
                    # compound = pcp.Compound.from_cid(int(row["CID"]))
                    compound = pcp.get_compounds(int(row["CID"]), "cid")
                    if row["new_name"] == "":
                        df.loc[idx, "new_name"] = compound[0].synonyms[0]
                except:
                    pass

            else:
                try:
                    compound = pcp.get_compounds(identifier=smi, namespace="smiles")
                    # if return more than one CID, will double check if the first one is good,
                    # else, write a comment
                    if len(compound) >= 2:
                        df.loc[idx, "comments"] = row["comment"] + \
                                                  "More than one CID available for this compound."
                    if row["new_name"] == "":
                        df.loc[idx, "new_name"] = compound[0].synonyms[0]
                    # add cid information
                    df.loc[idx, "CID"] = compound[0].cid
                except:
                    pass

            print(idx + 2)

    df.to_excel(out_excel, index=None)

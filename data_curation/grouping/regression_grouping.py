import pandas as pd
import numpy as np
from rdkit import Chem
from scipy import stats
import pubchempy as pcp

df = pd.read_excel("../2_bbb_all_complete_CID_out_smiles_fixed_updated.xlsx")
df = df[~df["logBB"].isna()]


df["logBB"] = df["logBB"].astype(float)
# remove molecules with logBB <= -9
df = df[df["logBB"] > -9]

# a dictionary to host inchi keys and isomeric smiles
for idx, row in df.iterrows():
    mol = Chem.MolFromSmiles(row["smiles_fixed_rdkit"])
    df.loc[idx, "Inchi"] = Chem.inchi.MolToInchi(mol)

df.to_excel("regression_inchi.xlsx", index=None, engine="openpyxl")
df = pd.read_excel("regression_inchi.xlsx")

# generate a dictionary to host all the inchi and isomeric smiles (or canonical smiles if isomeric smiles is not avaliable)


def append_value(dict_obj, key, value):
    if key in dict_obj:
        if not isinstance(dict_obj[key], list):
            dict_obj[key] = [dict_obj[key]]
        dict_obj[key].append(value)
    else:

        dict_obj[key] = value
    return dict_obj

inchi_smi_dict = {inchi:[] for inchi in df["Inchi"].to_list()}

for idx, row in df.iterrows():
    inchi_smi_dict = append_value(inchi_smi_dict, row["Inchi"], row["smiles_fixed_rdkit"])

# exam how inchi has more than one isomeric smiles
counter = 0
for key, value in inchi_smi_dict.items():
    if len(value) >= 2:
        counter += 1
print(counter)

# use non-redundant isomeric smiles for inchi_smi_dict
# manually inspect inchies with more than one non-redundant smiles
inchi_smi_dict = {inchi: set(smi) for inchi, smi in inchi_smi_dict.items()}
counter = 0
for key, value in inchi_smi_dict.items():
    if len(value) >= 2:
        print(key, value)
# the same inchi may have more than one inchi values, 12 in total
# but they are just resonance structure, so use inchi as an identifier

###########################################################################
df = pd.read_excel("regression_inchi.xlsx")

# smiles fixing with 02_clean_smiles_chembl_way_20210214.py


#########################################################################
df_unique = df.drop_duplicates(subset="Inchi", keep="first").reset_index(drop=True)
# df_duplicated = df.drop_duplicates(subset="Inchi", keep=False).reset_index(drop=True)

# df_unique["logBB"] = [[] for _ in np.arange(df_unique.shape[0])]
df_unique["logBB"] = ""
df_unique["compound_name"] = ""
df_unique["CID"] = ""
df_unique["new_name"] = ""
df_unique["iupac_name"] = ""
df_unique["reference"] = ""
df_unique["NO."] = ""

df["logBB"] = df["logBB"].astype(float)
# append compound_name, CID, logBB, new_name, iupac_name to the df_unique
# for idx_unique, row_unique in df_unique.iterrows():
#     for idx, row in df.iterrows():
#         if row["Inchi"] == row_unique["Inchi"]:
#             # logBB
#             df_unique.loc[idx_unique, "logBB"] = df_unique.loc[idx_unique, "logBB"] + "|" + str(row["logBB"])
#             # compound_name
#             df_unique.loc[idx_unique, "compound_name"] = df_unique.loc[idx_unique, "compound_name"] + "|" + str(row["compound_name"])
#             # CID
#             df_unique.loc[idx_unique, "CID"] = df_unique.loc[idx_unique, "CID"] + "|" + str(row["CID"])
#             # new_name
#             df_unique.loc[idx_unique, "new_name"] = df_unique.loc[idx_unique, "new_name"] + "|" + str(row["new_name"])
#             # iupac_name
#             df_unique.loc[idx_unique, "iupac_name"] = df_unique.loc[idx_unique, "iupac_name"] + "|" + str(row["iupac_name"])

# df_unique.to_excel("tmp.xlsx", index=None, engine="openpyxl")

# a more efficient way
for idx_unique, row_unique in df_unique.iterrows():
    inchi_unique = row_unique["Inchi"]
    df_inchi_matching = df[df["Inchi"] == inchi_unique].reset_index(drop=True)
    for _, row_matching in df_inchi_matching.iterrows():
        # logBB
        # df_unique.loc[idx_unique, "logBB"] = df_unique.loc[idx_unique, "logBB"] + str(row_matching["logBB"]) + "|"
        df_unique.loc[idx_unique, "logBB"] = df_unique.loc[idx_unique, "logBB"] + str(round(row_matching["logBB"], 2)) + "|"
        # compound_name
        df_unique.loc[idx_unique, "compound_name"] = df_unique.loc[idx_unique, "compound_name"] + str(row_matching["compound_name"]) + "|"
        # CID
        df_unique.loc[idx_unique, "CID"] = df_unique.loc[idx_unique, "CID"] + str(row_matching["CID"]) + "|"
        # new_name
        df_unique.loc[idx_unique, "new_name"] = df_unique.loc[idx_unique, "new_name"] + str(row_matching["new_name"]) + "|"
        # iupac_name
        df_unique.loc[idx_unique, "iupac_name"] = df_unique.loc[idx_unique, "iupac_name"] + str(row_matching["iupac_name"]) + "|"
        # reference
        df_unique.loc[idx_unique, "reference"] = df_unique.loc[idx_unique, "reference"] + str(row_matching["reference"]) + "|"
        # original NO.
        df_unique.loc[idx_unique, "NO."] = df_unique.loc[idx_unique, "NO."] + str(row_matching["NO."]) + "|"


df_unique.to_excel("regression_logBB_combined.xlsx", index=None, engine="openpyxl")

##################################################
# preprocess logBB data
from copy import deepcopy


df = pd.read_excel("regression_logBB_combined.xlsx")
# df_bak = deepcopy(df)

# filter molecules with max(logBB) – min(logBB) > 1
counter = 0
for idx, row in df.iterrows():
    logBB_values = [float(logBB) for logBB in row["logBB"].strip("|").split("|")]
    if max(logBB_values) - min(logBB_values) > 1:
        counter += 1
        df.loc[idx, "logBB"] = np.nan
df = df.dropna(subset=["logBB"]).reset_index(drop=True)

df["std"] = np.nan
df["group"] = ""

for idx, row in df.iterrows():
    # round logBB values to two decimal points as this is the most data hold for
    logBB_values = [logBB for logBB in row["logBB"].strip("|").split("|")]
    # find the minimum decimal places
    decimal_places = min([logBB[::-1].find('.') for logBB in logBB_values])
    logBB_values = [round(float(logBB), decimal_places) for logBB in logBB_values]
    # set logBB values if there is only one
    if len(logBB_values) == 1:
        df.loc[idx, "logBB"] = logBB_values[0]
        df.loc[idx, "group"] = "A"
        df.loc[idx, "std"] = 0
    else:
        mean_logBB = np.multiply(np.ones(len(logBB_values)),
                                 np.average(logBB_values))
        mean_logBB = np.around(mean_logBB, decimals=decimal_places)
        # set logBB values if all the values are the same or within 5% difference
        if np.allclose(np.array(logBB_values), mean_logBB, atol=0, rtol=0.05):
            df.loc[idx, "logBB"] = mean_logBB[0]
            df.loc[idx, "group"] = "B"
            df.loc[idx, "std"] = np.std(logBB_values)
        else:
            # if less than 3 values, use average value
            if len(logBB_values) < 3:
                df.loc[idx, "logBB"] = mean_logBB[0]
                df.loc[idx, "group"] = "C"
                df.loc[idx, "std"] = np.std(logBB_values)
            # if more than 3 values, use mode
            else:
                # not using stats.mode() because it can not handel the suitation when two mode values are avaliable
                # stats.mode(logBB_values)[0]
                values, counts = np.unique(logBB_values, return_counts=True)
                sorted_idx = np.argsort(counts)[::-1]
                values_sorted = values[sorted_idx]
                counts_sorted = counts[sorted_idx]
                # when there is only one number of maximum counts
                if counts_sorted[0] > counts_sorted[1]:
                    df.loc[idx, "logBB"] = values_sorted[0]
                    df.loc[idx, "group"] = "D"
                    df.loc[idx, "std"] = np.std(logBB_values)
                # when there are more than one maximum counts, they are equal
                else:
                    # more than 3 unique values
                    if len(values_sorted) >= 3:
                        # when there are two mode numbers
                        # counts_sorted[0] == counts_sorted[1] is a fact in such a condition as it
                        # is sorted
                        # the first 3 counts are the same
                        if counts_sorted[1] == counts_sorted[2]:
                            df.loc[idx, "logBB"] = sum(values_sorted[:3]) / 3
                            df.loc[idx, "group"] = "dropped_E"
                            df.loc[idx, "std"] = np.std(logBB_values)
                        # the first 2 counts are the same
                        else:
                            df.loc[idx, "logBB"] = sum(values_sorted[:2]) / 2
                            df.loc[idx, "group"] = "dropped_F"
                            df.loc[idx, "std"] = np.std(logBB_values)
                        # as counts_sorted is in descening order, counts_sorted[0] will not be less than counts_sorted[1]
                        # counts_sorted[0] == counts_sorted[1] and counts_sorted[0] == counts_sorted[2]
                    # when there are two unique count values
                    else:
                        # these two unique values are the same
                        if counts_sorted[0] == counts_sorted[1]:
                            df.loc[idx, "logBB"] = mean_logBB[0]
                            df.loc[idx, "group"] = "dropped_G"
                            df.loc[idx, "std"] = np.std(logBB_values)
                        # the first one is greater than the second one
                        else:
                            df.loc[idx, "logBB"] = values_sorted[0]
                            df.loc[idx, "group"] = "dropped_H"
                            df.loc[idx, "std"] = np.std(logBB_values)

#iupac name
for idx, row in df.iterrows():
    iupac_names = [name.lower() for name in row["iupac_name"].strip("|").split("|")
                   if name != "nan" if not name.isdigit() if len(name) != 1]
    if len(iupac_names) >= 1:
        df.loc[idx, "iupac_name"] = iupac_names[0].lstrip()
    else:
        df.loc[idx, "iupac_name"] = ""

# deal with compound_name, new_name
df["new_compound_name"] = ""
for idx, row in df.iterrows():
    # new_compound_name
    compound_names = [name.lower() for name in row["compound_name"].strip("|").split("|")
                      if name != "nan" if not name.isdigit() if len(name) != 1]
    new_names = [name.lower() for name in row["new_name"].strip("|").split("|")
                 if name != "nan" if not name.isdigit() if len(name) != 1]
    # these names found in pubchem come first
    names = list(set(new_names + compound_names))
    # when compound_names list is not empty
    if names != []:
        df.loc[idx, "new_compound_name"] = names[0].lstrip()
    else:
        df.loc[idx, "new_compound_name"] = row["iupac_name"]

# deal with CID
# for idx, row in df.iterrows():
#     cids = list(set([int(float(cid)) for cid in row["CID"].strip("|").split("|") if cid != "nan"]))
#     if len(cids) != 0:
#         df.loc[idx, "CID"] = cids[0]
#     else:
#         df.loc[idx, "CID"] = ""


# deal with smiles and CID
# df["smiles_fixed_rdkit"] = df["smiles_fixed_rdkit"].astype(str)
# df["CID"] = df["CID"].astype(str)

# for idx, row in df.iterrows():
#     # smiles_list = [smi.lower() for smi in row["smiles_fixed_rdkit"].strip("|").split("|")
#     #                if smi != "nan" if not smi.isdigit() if len(smi) != 1]

#     smiles_list = [smi for smi in row["smiles_fixed_rdkit"].strip("|").split("|")
#                    if smi != "nan" if not smi.isdigit()]
#     smiles_list = list(set(smiles_list))

#     cids = list(set([int(float(cid)) for cid in row["CID"].strip("|").split("|") if cid != "nan"]))

#     if len(smiles_list) >= 1:
#         # df.loc[idx, "smiles_fixed_rdkit"] = smiles_list[0].lstrip()

#         # get new CID from the smiles if CID is none
#         # else: use old CID
#         if len(cids) == 0:
#             ## try to get CID until using up the smiles

#             # flag to indicate if we found new CID and smiles
#             flag = False
#             for smi in smiles_list:
#                 try:
#                     # because can get an error with
#                     # O=[SH](O)(c1ccc2cc[nH]c2c1)N1CCCC1CCN1CCC(Oc2cccc(Cl)c2)CC1
#                     compound = pcp.get_compounds(identifier=smi, namespace="smiles")
#                     cid_new = compound[0].cid
#                     if cid_new is not None:
#                         flag = True
#                         break
#                 except:
#                     print("error found when searching pubchem")

#             if flag is True:
#                 df.loc[idx, "smiles_fixed_rdkit"] = smi
#                 df.loc[idx, "CID"] = cid_new
#             else:
#                 df.loc[idx, "smiles_fixed_rdkit"] = smiles_list[0]
#                 df.loc[idx, "CID"] = ""
#         else:
#             # use old CIDs
#             df.loc[idx, "smiles_fixed_rdkit"] = smiles_list[0]
#             if len(cids) >= 1:
#                 df.loc[idx, "CID"] = cids[0]
#             else:
#                 df.loc[idx, "CID"] = ""
###########################################################
df["CID"] = df["CID"].fillna("")
df["CID"] = df["CID"].astype(str)

# deal with CID
df["CID"] = df["CID"].astype(str)
for idx, row in df.iterrows():
    # no need to deal with CID for regression data again
    if pd.isnull(row["logBB"]):
        cids = list(set([int(float(cid)) for cid in row["CID"].strip("|").split("|") if cid != "nan"]))
        if len(cids) != 0:
            df.loc[idx, "CID"] = cids[0]
        else:
            df.loc[idx, "CID"] = ""


# deal with SMILES
df["smiles_fixed_rdkit"] = df["smiles_fixed_rdkit"].astype(str)
for idx, row in df.iterrows():
    smi_strings = list(set([smi for smi in row["smiles_fixed_rdkit"].strip("|").split("|") if smi != "nan"]))
    if len(cids) != 0:
        df.loc[idx, "smiles_fixed_rdkit"] = smi_strings[0]
    else:
        df.loc[idx, "smiles_fixed_rdkit"] = ""


df = df.sort_values(by=["group", "logBB"])
df.to_excel("regression_clean_done.xlsx", index=None, engine="openpyxl")
# clean the data manually

import pandas as pd
import numpy as np
from rdkit import Chem
from scipy import stats
import pubchempy as pcp

df = pd.read_excel("../2_bbb_all_complete_CID_out_smiles_fixed_updated.xlsx")

df["logBB"] = df["logBB"].astype(float)
# filter out those values with logBB <= -9
df = df[(df["logBB"] > -9) | (df["logBB"].isna())]

# a dictionary to host inchi keys and isomeric smiles
for idx, row in df.iterrows():
    mol = Chem.MolFromSmiles(row["smiles_fixed_rdkit"])
    df.loc[idx, "Inchi"] = Chem.inchi.MolToInchi(mol)

df.to_excel("classification_inchi.xlsx", index=None, engine="openpyxl")
df = pd.read_excel("classification_inchi.xlsx")

# generate a dictionary to host all the inchi and isomeric smiles (or canonical smiles if isomeric smiles is not avaliable)


def append_value(dict_obj, key, value):
    if key in dict_obj:
        if not isinstance(dict_obj[key], list):
            dict_obj[key] = [dict_obj[key]]
        dict_obj[key].append(value)
    else:

        dict_obj[key] = value
    return dict_obj

# only exam those without logBB values
# because those with logBB have been examined in regression
df_tmp = df[df["logBB"].isna()]
inchi_smi_dict = {inchi:[] for inchi in df_tmp["Inchi"].to_list()}

for idx, row in df_tmp.iterrows():
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

###################################################################
df_original = pd.read_excel("classification_inchi.xlsx")

# build a new dataframe with regression type and classification type
# df_classification_only = df_original.loc[df_original["logBB"].isna()]

df_classification_only = df_original[df_original["logBB"].isna()]
df_regression_only = pd.read_excel("../regression/regression_clean_done.xlsx")
df_regression_only = df_regression_only[~df_regression_only["group"].str.contains("dropped")]

# add threshold values for df_classification_only
threshold_dict = {
        "R1": "-1",
        "R2": "(-2,-1)",
        "R7": "-1",
        "R8": "0",
        "R13": "-1",
        "R16": "-1",
        "R18": "0.1",
        "R19": "-1",
        "R25": "0",
        "R28": "-1"
        }

# drop threshold column
# df_classification_only.drop("threshold", axis=1, inplace=True)
# df_classification_only["threshold"] = df_classification_only.apply(get_threshold, axis=1)
for idx, row in df_classification_only.iterrows():
    df_classification_only.loc[idx, "threshold"] = threshold_dict.get(row["reference"])

# sort df_classification according the threshold value
# -1, 0, 0.1, (-2,1)
# but the fact is that only -1 exists for categorical data
# other threshold values are in the regression data type

# add a new columns to set as a sorting column
# df["reference"] = df["reference"].astype(str)
df_classification_only["sort_ref"] = df_classification_only["reference"].apply(lambda x: int(x.lstrip("R")))
df_classification_only = df_classification_only.sort_values(by=["threshold", "sort_ref"], ascending=[True, True], inplace=False)
df_classification_only.drop("sort_ref", axis=1, inplace=True)
df_classification_only.to_excel("df_classification_only.xlsx", index=None, engine="openpyxl")


# contact regression and classification together
df_all = pd.concat([df_regression_only, df_classification_only], axis=0)
df_all.to_excel("df_classification_reformatted.xlsx", index=None, engine="openpyxl")



df_all = pd.read_excel("df_classification_reformatted.xlsx")
df_all = df_all.sort_values(by=["logBB"]).reset_index(drop=True)
# unique records by inchi values
df_unique = df_all.drop_duplicates(subset="Inchi", keep="first").reset_index(drop=True)

# group A: records with logBB values
# df_unique["logBB"] = [[] for _ in np.arange(df_unique.shape[0])]
df_unique["smiles_fixed_rdkit"] = ""
df_unique["BBB+/BBB-"] = ""
df_unique["compound_name"] = ""
df_unique["CID"] = ""
df_unique["new_name"] = ""
df_unique["iupac_name"] = ""
df_unique["reference"] = ""
df_unique["NO."] = ""



# df_unique["logBB"] = df_unique["logBB"].astype(float)

# deal with data without logBB values
for idx_unique, row_unique in df_unique.iterrows():
    inchi_unique = row_unique["Inchi"]
    df_inchi_matching = df_classification_only[df_classification_only["Inchi"] == inchi_unique].reset_index(drop=True)
    for _, row_matching in df_inchi_matching.iterrows():
        # SMILES
        df_unique.loc[idx_unique, "smiles_fixed_rdkit"] = df_unique.loc[idx_unique, "smiles_fixed_rdkit"] + str(row_matching["smiles_fixed_rdkit"]) + "|"
        # BBB+/BBB-
        df_unique.loc[idx_unique, "BBB+/BBB-"] = df_unique.loc[idx_unique, "BBB+/BBB-"] + str(row_matching["BBB+/BBB-"]) + "|"
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


# df_combined = pd.concat([df_regression_only, df_unique[df_unique["logBB"].isna()]])
# # df_combined.to_excel("tmp_combined.xlsx", index=None, engine="openpyxl")
# df_combined.to_excel("classification_logBB_combined.xlsx", index=None, engine="openpyxl")
df_unique.to_excel("classification_logBB_combined.xlsx", index=None, engine="openpyxl")

##################################################
# preprocess BBB+/BBB- data
from copy import deepcopy


df = pd.read_excel("classification_logBB_combined.xlsx")
# df_bak = deepcopy(df)

# filter molecules with max(logBB) – min(logBB) > 1
# counter = 0
# for idx, row in df.iterrows():
#     logBB_values = [float(logBB) for logBB in row["logBB"].strip("|").split("|")]
#     if max(logBB_values) - min(logBB_values) > 1:
#         counter += 1
#         df.loc[idx, "logBB"] = np.nan
# df = df.dropna(subset=["logBB"]).reset_index(drop=True)

df["group"] = ""

threshold_logBB = -1

for idx, row in df.iterrows():
    # group A: those with logBB values
    if not pd.isnull(row["logBB"]):
        # assign logBB values
        if row["logBB"] >= threshold_logBB:
            df.loc[idx, "BBB+/BBB-"] = "BBB+"
        else:
            df.loc[idx, "BBB+/BBB-"] = "BBB-"
        df.loc[idx, "group"] = "A"
    else:
        categories = [cat for cat in row["BBB+/BBB-"].strip("|").split("|")]
        categories_set = list(set(categories))
        # threshold = -1
        if row["threshold"] == threshold_logBB:
            # group B: threshold=-1 and consistent results
            if len(categories_set) == 1:
                df.loc[idx, "BBB+/BBB-"] = categories[0]
                df.loc[idx, "group"] = "B"
            # group C: threshold=-1 and inconsistent results
            # 203 records
            else:
                # old group C -> dropped
                df.loc[idx, "group"] = "dropped_group1"
                df.loc[idx, "BBB+/BBB-"] = ""
        # those without threshold values
        else:
            # group D: threshold not given, consistent results
            if len(categories_set) == 1:
                df.loc[idx, "BBB+/BBB-"] = categories_set[0]
                # old group D -> C
                df.loc[idx, "group"] = "C"
            else:
                values, counts = np.unique(categories, return_counts=True)
                sorted_idx = np.argsort(counts)[::-1]
                values_sorted = values[sorted_idx]
                counts_sorted = counts[sorted_idx]
                # group E: threshold not given, inconsistent results, but voting works
                if counts_sorted[0] > counts_sorted[1]:
                    df.loc[idx, "BBB+/BBB-"] = values_sorted[0]
                    # old group E -> D
                    df.loc[idx, "group"] = "D"
                # group F: threshold not given, inconsistent results, voting doesn't work
                else:
                    # old group F -> dropped
                    df.loc[idx, "group"] = "dropped_group2"
                    df.loc[idx, "BBB+/BBB-"] = ""


df.to_excel("tmp.xlsx", index=None, engine="openpyxl")


df = pd.read_excel("tmp.xlsx")
# df = pd.read_excel("tmp_test.xlsx")
df = df[~df["group"].str.contains("dropped")]

# iupac name
df["iupac_name"] = df["iupac_name"].astype(str)
for idx, row in df.iterrows():
    iupac_names = [name.lower() for name in row["iupac_name"].strip("|").split("|")
                   if name != "nan" if not name.isdigit() if len(name) != 1]
    if len(iupac_names) >= 1:
        df.loc[idx, "iupac_name"] = iupac_names[0].lstrip()
    else:
        df.loc[idx, "iupac_name"] = ""

# deal with compound_name, new_name
df["new_compound_name"] = ""
df["compound_name"] = df["compound_name"].astype(str)
df["new_name"] = df["new_name"].astype(str)

for idx, row in df.iterrows():
    # new_compound_name
    # compound_names = [name.lower() for name in row["compound_name"].strip("|").split("|")
    #                   if name != "nan" if not name.isdigit() if len(name) != 1]
    # new_names = [name.lower() for name in row["new_name"].strip("|").split("|")
    #              if name != "nan" if not name.isdigit() if len(name) != 1]

    compound_names = [name.lower() for name in row["compound_name"].strip("|").split("|")
                      if name != "nan" if not name.isdigit()]
    new_names = [name.lower() for name in row["new_name"].strip("|").split("|")
                 if name != "nan" if not name.isdigit()]

    # these names found in pubchem come first
    names = list(set(new_names + compound_names))
    # when compound_names list is not empty
    if names != []:
        df.loc[idx, "new_compound_name"] = names[0].lstrip()
    else:
        df.loc[idx, "new_compound_name"] = row["iupac_name"]


# deal with smiles and CID
# df["smiles_fixed_rdkit"] = df["smiles_fixed_rdkit"].astype(str)

df["CID"] = df["CID"].fillna("")
df["CID"] = df["CID"].astype(str)

################################################
# # old way of doing
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
############################################################

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
df.to_excel("classification_clean_done.xlsx", index=None, engine="openpyxl")


################################################################

# clean the data manually

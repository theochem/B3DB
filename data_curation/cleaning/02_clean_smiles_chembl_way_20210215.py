#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modified on Sun Feb 9 2021. Molecule cleaning for BBB database.

@author: Fanwang Meng @ Ayers Lab
Usage:
"""

import argparse

import pandas as pd
from chembl_structure_pipeline.exclude_flag import exclude_flag
from chembl_structure_pipeline.standardizer import (
    get_parent_mol, update_mol_valences, kekulize_mol, flatten_tartrate_mol,
    normalize_mol, uncharge_mol, cleanup_drawing_mol)
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolops


metal_atoms = {3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,
               29, 30, 31, 32, 34, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46,
               47, 48, 49, 50, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65,
               66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80,
               81, 82, 83, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98,
               99, 100, 101, 102, 103}
heavy_atoms = {21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34,
               37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
               50, 51, 52, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64,
               65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77,
               78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90,
               91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102,
               103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113,
               114, 115, 116, 117, 118}

unwanted_atom_dict = {
    "metal_atoms": metal_atoms,
    "heavy_atoms": heavy_atoms,
}


class CleanMoleculesFromDataFrame:
    def __init__(self,
                 input_fname="2_bbb_all_complete_CID_out_smiles.xlsx",
                 output_fname="2_bbb_all_complete_CID_out_smiles_fixed.xlsx",
                 unwanted_atoms="heavy_atoms",
                 ):
        """
        Clean molecules from a single SDF file with multiple molecules.

        Notes
        -----
        1. Now we only support SDF file format only. For MOL2 format, one may refer to
        https://chem-workflows.com/articles/2020/03/23/building-a-multi-molecule-mol2-reader-for
        -rdkit-v2/, which remains to be tested.
        2. When using SaltRemover from rdkit, all the hydrogen should not be imported, but should be
        ignored. Otherwise, it will not work.

        """
        self.input_fname = input_fname
        self.output_fname = output_fname
        self.unwanted_atoms = unwanted_atoms

        self.df = pd.read_excel(self.input_fname)

    @staticmethod
    def skip_problematic_molecules(df,
                                   second_check=False,
                                   ignore_none=True,
                                   ignore_empty_mol=False):
        """Skip problematic molecules."""
        if second_check:
            smi_label = "smiles_fixed_rdkit"
        else:
            smi_label = "smiles_result"

        for idx, row in df.iterrows():
            try:
                mol = Chem.MolFromSmiles(row[smi_label])
                if ignore_none is True and ignore_empty_mol is False:
                    if mol is not None:
                        df.loc[idx, "smiles_fixed_rdkit"] = Chem.MolToSmiles(mol, canonical=True)
                    else:
                        df.loc[idx, "smiles_fixed_rdkit"] = ""

                if ignore_none is False and ignore_empty_mol is True:
                    if mol.GetNumAtoms() != 0:
                        df.loc[idx, "smiles_fixed_rdkit"] = Chem.MolToSmiles(mol, canonical=True)
                    else:
                        df.loc[idx, "smiles_fixed_rdkit"] = ""

                if ignore_none is True and ignore_empty_mol is True:
                    if mol is not None and mol.GetNumAtoms() != 0:
                        df.loc[idx, "smiles_fixed_rdkit"] = Chem.MolToSmiles(mol, canonical=True)
                    else:
                        df.loc[idx, "smiles_fixed_rdkit"] = ""

            except:
                df.loc[idx, "smiles_fixed_rdkit"] = ""
                print(row["NO."])


        # df = df.dropna(axis=0, subset=["smiles_fixed_rdkit"], inplace=False)
        df = df[df["smiles_fixed_rdkit"] != ""]
        df = df.reset_index(drop=True)
        return df

    @staticmethod
    def filter_restricted_atoms(df, unwanted_atoms="heavy_atoms"):
        """Filter by unwanted atom set."""
        unwanted_atom_set = unwanted_atom_dict.get(unwanted_atoms)

        metal_count = 0
        for idx, row in df.iterrows():
            if row["smiles_fixed_rdkit"] != "":
                try:
                    mol = Chem.MolFromSmiles(row["smiles_fixed_rdkit"])
                    atom_list = [atom.GetAtomicNum() for atom in mol.GetAtoms()]

                    if not set(atom_list).isdisjoint(unwanted_atom_set):
                        df.loc[idx, "smiles_fixed_rdkit"] = ""
                        metal_count += 1
                except AttributeError:
                    continue

        df = df[df["smiles_fixed_rdkit"] != ""]
        df = df.reset_index(drop=True)
        return df

    #######################################################
    # # get_parent_mol() already implemented this functionality
    # @staticmethod
    # def keep_largest_fragment(df):
    #     """Only keep the largest molecular fragment"""

    #     for idx, row in df.iterrows():
    #         if row["smiles_fixed_rdkit"] != "":
    #             mol = Chem.MolFromSmiles(row["smiles_fixed_rdkit"])
    #             mols = rdmolops.GetMolFrags(mol, asMols=True)
    #             mol_new = max(mols, default=mol, key=lambda m: m.GetNumAtoms())
    #             df.loc[idx, "smiles_fixed_rdkit"] = Chem.MolToSmiles(mol_new,
    #                                                                  canonical=True)

    #     return df

    @staticmethod
    def handle_solvents_salts(df,
                              neutralize=False,
                              check_exclusion=True,
                              verbose=False):
        """"""
        # get_parent_mol() function
        # can lead to an error of converting smiles
        # if input smiles to build a molecule is
        # c1(oc([CH2+]SCCN\C(=C\[N+](=O)[O-])NC)cc1)CN(C)C
        # output molecule smiles would be
        # CN/C(=C\\[N+](=O)[O-])NCCS[CH3]c1ccc(CN(C)C)o1
        # however, loading CN/C(=C\\[N+](=O)[O-])NCCS[CH3]c1ccc(CN(C)C)o1 cannot build a molecule
        # df['column'] = df['column'].fillna(value)
        df["smiles_fixed_rdkit"] = df["smiles_fixed_rdkit"].fillna("")

        for idx, row in df.iterrows():
            if row["smiles_fixed_rdkit"] != "":
                mol = Chem.MolFromSmiles(row["smiles_fixed_rdkit"])
                mol_new, _ = get_parent_mol(mol,
                                            neutralize=neutralize,
                                            check_exclusion=check_exclusion,
                                            verbose=verbose)

                df.loc[idx, "smiles_fixed_rdkit"] = Chem.MolToSmiles(mol_new,
                                                                     canonical=True)

        df = df[df["smiles_fixed_rdkit"] != ""]
        df = df.reset_index(drop=True)
        return df

    @staticmethod
    def standardize_molecules(df, check_exclusion=True):
        """Standardize molecules with ChEMBL Structure Pipeline."""
        for idx, row in df.iterrows():
            if row["smiles_fixed_rdkit"] != "":
                mol = Chem.MolFromSmiles(row["smiles_fixed_rdkit"])

                if check_exclusion:
                    # Rules to exclude structures.
                    # Metallic or non metallic with more than 7 boron atoms will be excluded
                    # due to problems when depicting borane compounds.
                    exclude = exclude_flag(mol, includeRDKitSanitization=False)
                else:
                    exclude = False

                if not exclude:
                    mol = update_mol_valences(mol)
                    # removes all SubstanceGroups from a molecule (if any)
                    # mol = remove_sgroups_from_mol(mol)

                    mol = kekulize_mol(mol)

                    # mol = remove_hs_from_mol(mol)

                    # https://github.com/chembl/ChEMBL_Structure_Pipeline/blob
                    # /badad4f0432e7139adaa2ddac9e1c081be1e90fd/chembl_structure_pipeline
                    # /standardizer.py#L38-62
                    mol = normalize_mol(mol)
                    mol = uncharge_mol(mol)
                    # mol = flatten_tartrate_mol(mol)
                    # cleanup_drawing_mol() only works for 2D molecules
                    # mol = cleanup_drawing_mol(mol)

                df.loc[idx, "smiles_fixed_rdkit"] = Chem.MolToSmiles(mol, canonical=True)

        df = df[df["smiles_fixed_rdkit"] != ""]
        df = df.reset_index(drop=True)
        return df

    @staticmethod
    def neutralize(df):
        # Neutralizing molecules
        for idx, row in df.iterrows():
            if row["smiles_fixed_rdkit"] != "":
                smi_new, _ = neutralize_charges(smiles=row["smiles_fixed_rdkit"], reactions=None)
                df.loc[idx, "smiles_fixed_rdkit"] = smi_new

        return df

    def __call__(self):
        # filter molecules by if readable with RdKit and if unwanted atoms in molecules

        df = self.df
        df["smiles_fixed_rdkit"] = ""
        df = self.skip_problematic_molecules(df,
                                             ignore_none=True,
                                             second_check=False,
                                             ignore_empty_mol=False)

        # df = self.filter_restricted_atoms(df,
        #                                   self.unwanted_atoms)

        # keep the largest fragment
        # df = self.keep_largest_fragment(df)

        # get parent molecule by removing isotopes, solvents and salts
        # molecule neutralization is disabled for now as it can lead to an error
        df = self.handle_solvents_salts(df,
                                        neutralize=False,
                                        check_exclusion=True,
                                        verbose=False)

        df = self.neutralize(df)

        df = self.filter_restricted_atoms(df,
                                          self.unwanted_atoms)

        # standardize molecules including uncharge molecules
        df = self.standardize_molecules(df, check_exclusion=True)

        # skip problematic molecules
        df = self.skip_problematic_molecules(df,
                                             ignore_none=True,
                                             second_check=True,
                                             ignore_empty_mol=True)
        # sort result with order: logBB, CID, ...
        df = df.sort_values(by=["logBB", "smiles_fixed_rdkit", "CID"])
        # save result
        df.to_excel(self.output_fname, index=None, engine="openpyxl")


_reactions = None


def neutralize_charges(smiles, reactions=None):
    """ contribution from Hans de Winter """
    global _reactions
    if reactions is None:
        if _reactions is None:
            _reactions = _InitialiseNeutralisationReactions()
        reactions = _reactions
    mol = Chem.MolFromSmiles(smiles)
    replaced = False
    for i, (reactant, product) in enumerate(reactions):
        while mol.HasSubstructMatch(reactant):
            replaced = True
            rms = AllChem.ReplaceSubstructs(mol, reactant, product)
            mol = rms[0]
    if replaced:
        return Chem.MolToSmiles(mol, True), True
    else:
        return smiles, False


def _InitialiseNeutralisationReactions():
    """Contribution from Hans de Winter."""
    patts = (
        # Imidazoles
        ('[n+;H]', 'n'),
        # Amines
        ('[N+;!H0]', 'N'),
        # Carboxylic acids and alcohols
        ('[$([O-]);!$([O-][#7])]', 'O'),
        # Thiols
        ('[S-;X1]', 'S'),
        # Sulfonamides
        ('[$([N-;X2]S(=O)=O)]', 'N'),
        # Enamines
        ('[$([N-;X2][C,N]=C)]', 'N'),
        # Tetrazoles
        ('[n-]', '[nH]'),
        # Sulfoxides
        ('[$([S-]=O)]', 'S'),
        # Amides
        ('[$([N-]C=O)]', 'N'),
        )
    return [(Chem.MolFromSmarts(x), Chem.MolFromSmiles(y, False)) for x, y in patts]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-input_fname", type=str,
                        default="2_bbb_all_complete_CID_out_smiles.xlsx",
                        help="Input excel file.")
    parser.add_argument("-output_fname", type=str,
                        default="2_bbb_all_complete_CID_out_smiles_fixed.xlsx",
                        help="File to output results.")
    parser.add_argument("-unwanted_atoms", type=str,
                        default="heavy_atoms",
                        help="Unwanted atom set used to filter molecules.")
    args = parser.parse_args()

    cleaning_molecules = CleanMoleculesFromDataFrame(
        input_fname=args.input_fname,
        output_fname=args.output_fname,
        unwanted_atoms=args.unwanted_atoms,
    )
    cleaning_molecules()

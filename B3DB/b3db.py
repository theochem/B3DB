# The B3DB provides a rich molecule database for
# Blood-Brain Barrier (BBB) permeability.
#
# Copyright (C) 2021-2025 The QC-Devs Community
# This file is part of B3DB.
#
# B3DB is dedicated to the public domain under the terms of CC0 1.0 Universal.
#
# A full copy of the CC0 1.0 Universal license can be found at:
# https://creativecommons.org/publicdomain/zero/1.0/legalcode
#
# --

from pathlib import Path

import pandas as pd


def load_b3db_dataset():
    """Load the B3DB dataset."""
    data_dir = Path(__file__).parent
    # load regression dataset
    regression_data = pd.read_csv(data_dir / "B3DB_regression.tsv", sep="\t")
    # load extended regression dataset
    regression_data_extended = pd.read_csv(
        data_dir / "B3DB_regression_extended.tsv.gz", sep="\t", compression="gzip", low_memory=False
    )
    # load classification dataset
    classification_data = pd.read_csv(data_dir / "B3DB_classification.tsv", sep="\t")
    # load extended classification dataset
    classification_data_extended = pd.read_csv(
        data_dir / "B3DB_classification_extended.tsv.gz",
        sep="\t",
        compression="gzip",
        low_memory=False,
    )
    # load manually curated external data
    classification_external_data = pd.read_csv(
        data_dir / "B3DB_classification_external.tsv", sep="\t"
    )

    return {
        "B3DB_regression": regression_data,
        "B3DB_classification": classification_data,
        "B3DB_regression_extended": regression_data_extended,
        "B3DB_classification_extended": classification_data_extended,
        "B3DB_classification_external": classification_external_data,
    }


B3DB_DATA_DICT = load_b3db_dataset()

__all__ = [
    "B3DB_DATA_DICT",
]

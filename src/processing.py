from typing import List

from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw

import numpy as np


def extract_smiles_from_file(path_to_smi: str) -> tuple[List, List]:
    """
    Extract smiles from smi file

    Args:
        path_to_smi:

    Returns:

    """

    try:
        with open(path_to_smi, "r") as smi_file:
            smiles_name_list = smi_file.read().splitlines()
    except FileNotFoundError:
        return [], []

    names_list = []
    smiles_list = []
    for sm in smiles_name_list:
        try:
            smile, name = sm.split("#")
        except ValueError:
            smile = sm
            name = ""
        names_list.append(name)
        smiles_list.append(smile)

    return smiles_list, names_list


def mol_to_numpy(smiles):
    """
    Convert mol to numpy

    Args:
        smiles:

    Returns:

    """
    mol = Chem.MolFromSmiles(smiles)
    img = np.array(Draw.MolToImage(mol))
    return img

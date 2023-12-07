from rdkit.Chem import Descriptors
from rdkit import Chem

from .processing import extract_smiles_from_file


def generate_property(path_to_smi: str) -> dict:
    """
    Generate property from smiles

    Args:
        path_to_smi:

    Returns:

    """
    smiles_list, names_list = extract_smiles_from_file(path_to_smi)

    property_dict = {}
    for smiles, name in zip(smiles_list, names_list):
        mol = Chem.MolFromSmiles(smiles)

        mol_weight = Descriptors.MolWt(mol)
        mol_logP = Descriptors.MolLogP(mol)  # noqa
        mol_num_ha_acceptors = Descriptors.NumHAcceptors(mol)
        mol_num_hd_donors = Descriptors.NumHDonors(mol)

        property_dict[smiles] = {
            "name": name,
            "weight": mol_weight,
            "logP": mol_logP,
            "num_ha_acceptors": mol_num_ha_acceptors,
            "num_hd_donors": mol_num_hd_donors,
        }

    return property_dict

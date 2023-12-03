# %%
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from typing import Optional
import PIL
from PIL import Image


def generate(
    smiles_list,
    mols_per_row: Optional[int] = 3,
    sub_img_size: Optional[tuple] = (200, 200),
    add_smiles_name: Optional[bool] = False,
):
    mols = [Chem.MolFromSmiles(m) for m in smiles_list]

    # display
    if add_smiles_name:
        img = Draw.MolsToGridImage(
            mols,
            molsPerRow=mols_per_row,
            subImgSize=sub_img_size,
            legends=[murckoHash for murckoHash in smiles_list],
            useSVG=False,
            returnPNG=False,
        )
    else:
        img = Draw.MolsToGridImage(
            mols,
            molsPerRow=mols_per_row,
            subImgSize=sub_img_size,
            useSVG=False,
            returnPNG=False,
        )
    return img


def processing_smi(
    path_to_smi,
    mols_per_row: Optional[int] = 3,
    sub_img_size: Optional[tuple] = (200, 200),
    add_smiles_name: Optional[bool] = False,
):
    with open(path_to_smi, "r") as smi_file:
        smiles_list = smi_file.read().splitlines()

    img = generate(
        smiles_list,
        mols_per_row=mols_per_row,
        sub_img_size=sub_img_size,
        add_smiles_name=add_smiles_name,
    )

    return img


def save_img(img: PIL.PngImagePlugin.PngImageFile, path_to_save: str):
    img.save(path_to_save)

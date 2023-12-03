# %%
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from typing import Optional, List
import PIL
from PIL import Image


def generate(
    smiles_list: List[str],
    names_list: List[str],
    mols_per_row: Optional[int] = 3,
    sub_img_size: Optional[tuple] = (200, 200),
    add_smiles_name: Optional[str] = "Not Annotation",
):
    mols = [Chem.MolFromSmiles(m) for m in smiles_list]
    print(mols)
    print(mols_per_row)
    print(sub_img_size)
    print(smiles_list)
    # display
    if add_smiles_name == "Smiles Annotation":
        img = Draw.MolsToGridImage(
            mols,
            molsPerRow=mols_per_row,
            subImgSize=sub_img_size,
            legends=[murckoHash for murckoHash in smiles_list],
            useSVG=False,
            returnPNG=False,
        )
    elif add_smiles_name == "Name Annotation":
        img = Draw.MolsToGridImage(
            mols,
            molsPerRow=mols_per_row,
            subImgSize=sub_img_size,
            legends=[name for name in names_list],
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
    path_to_smi: str,
    path_to_save: str,
    add_smiles_name: str,
    mols_per_row: Optional[int] = 3,
    sub_img_size: Optional[tuple] = (200, 200),
):
    with open(path_to_smi, "r") as smi_file:
        smiles_name_list = smi_file.read().splitlines()

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

    print(names_list)

    img = generate(
        smiles_list,
        names_list,
        mols_per_row=mols_per_row,
        sub_img_size=sub_img_size,
        add_smiles_name=add_smiles_name,
    )

    save_img(img, path_to_save)


def save_img(img: PIL.PngImagePlugin.PngImageFile, path_to_save: str):
    img.save(path_to_save)

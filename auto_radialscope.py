import pandas as pd
from get_r_groups import get_r_df, r_groups_from_scaffold
from get_r_data import get_smiles_and_activity
from calculate_drawing_angle import get_bond_angle, create_r_group_map
from rdkit import Chem
from rdkit.Chem import PandasTools


default_settings = {'SMILESSTRING': '',  # use methyl groups for your rests, smiles strings can't use labels
                    'use_bw_atom_theme': False,  # draw all atoms black
                    'use_bold_font': True, # replace all fonts wherever possible with bold text
                    'white_cutoff': 80,  #make text of labels white if greater than this number
                    'scalefactor': 0.8  #scale the total plot in the end, for large molecule you will need to decrease this or make the viewbox bigger in a vector software such as Inkscape
                    }

def make_auto_radial_scope(enzyme_name, df, scaffold_smi, activity_col='conversion'):

    scaffold = Chem.MolFromSmiles(scaffold_smi)

    # 1. Filter df by enzyme name, remove nan activity and any duplicate products
    # (it may be necessary to do some filters before this, eg for a single paper)
    df = df[df['enzyme_name'] == enzyme_name]
    df = df.dropna(subset=[activity_col])
    df = df.drop_duplicates(subset=['product_1_smiles'], keep=False)

    # 2. get r groups in a df and merge
    r_groups = r_groups_from_scaffold(scaffold)
    r_group_df = get_r_df(df, scaffold)
    df = pd.merge(df, r_group_df, on='product_1_smiles')

    # 3. get smiles and activity for radial scope
    #   - r_smiles will be a dict like {'R1': ['CCO', 'Br', ...], }
    #   - r_conversions will be a dict like {'R1': [90, 54, ...] }
    r_smiles, r_activity = get_smiles_and_activity(df, r_groups, activity_col=activity_col)

    # 4. get radial_scope_template_from_scaffold
    r_group_map = create_r_group_map(scaffold)
    settings = default_settings
    settings['SMILESSTRING'] = scaffold_smi


if __name__ == "__main__":
    # Load the data and filter as in Mengyu's notebook
    df = pd.read_excel('AL data.xlsx')
    df = df.replace("Œ±-amino amination", "alpha-amino amination")
    df = df.replace("Œ≤-amino deamination", "beta-amino amination")
    df = df[df['enzyme_type'] == "PAL"]
    df = df.dropna(subset=['product_1_smiles'])
    PandasTools.AddMoleculeColumnToFrame(df, 'product_1_smiles', 'Product Mol')

    scaffold_smi = 'N[C@@H](CC1=C([*:1])C([*:2])=C([*:3])C=C1)C(O)=O'
    make_auto_radial_scope("EncP_E293M", df, scaffold_smi)












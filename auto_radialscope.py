import pandas as pd
from get_r_groups import get_r_df, r_groups_from_scaffold
from get_r_data import get_smiles_and_activity, combine_enzyme_dict
from calculate_drawing_angle import get_bond_angle, create_r_group_map, calculate_r_scope_angle
from rdkit import Chem
from rdkit.Chem import PandasTools
from radialscope import RadialScope as rs
import copy

default_settings = {'SMILESSTRING': '',  # this will be updated
                    'use_bw_atom_theme': False,  # draw all atoms black
                    'use_bold_font': True,  # replace all fonts wherever possible with bold text
                    'white_cutoff': 80,  # make text of labels white if greater than this number
                    'scalefactor': 0.8
                    # scale the total plot in the end, for large molecule you will need to decrease this or make the viewbox bigger in a vector software such as Inkscape
                    }

default_r_setup = {'rest_label': "",  # Label of the atom in the radial scope plot eg R$_1$
                    'no_wedges': 0,  # Number of wedges in the Radial scope Plot
                    'coverangle_wedges': 0,  # 15*num_wedges - Degrees of a full circle the Rscope should cover
                    'startangle': 0,  # Start angle of the Rscope
                    'CMAPINNER': "Blues",  # Colormap of the inner circle. Max value has full color, options see below
                    'CMAPOUTER': "Greens",  # Colormap of the outer circle. Max value has full color
                    'OUTERLABEL': "",  # eg EncP_E293M % conversion
                    'INNERLABEL': "",  # Label of the inner circle
                    'value_inner_circle': [], # conversions
                    'value_outer_circle': [],
                    'min_max_value': [(0, 0), (100, 100)],
                    # Define minimum and maximum for each colorbar [(inner_min, outer_min),(inner_max, outer_max)]
                    'rounding': True,  # ENABLE rounding if you want something like >99
                    'rounding_boundary': 99,  # cutoff for displaying >
                    'value_groups': [], # labels for smiles
                    # Labels for the outer circle, you can use math using $, any ~ will be interpreted as smiles string
                    'attach_atom_id': 0,
                }


def get_r_setup_dict(r_num, enzyme_data, enzyme_1, enzyme_2, scaffold, r_group_map, angle_per_wedge=15):
    r_dict = copy.deepcopy(default_r_setup)
    r_string = f"R{r_num}"
    num_wedges = len(enzyme_data[r_string])
    if num_wedges == 0:
        return None

    r_dict['rest_label'] = f"R$_{r_num}$"  # eg "R$_1$"
    r_dict['no_wedges'] = num_wedges
    r_dict['coverangle_wedges'] = angle_per_wedge * num_wedges
    bond_angle = get_bond_angle(scaffold, r_num)
    r_dict['startangle'] = calculate_r_scope_angle(bond_angle, num_wedges, angle_per_wedge)
    r_dict['INNERLABEL'] = f"{enzyme_1} % conversion"
    r_dict['OUTERLABEL'] = f"{enzyme_2} % conversion"

    r_dict['value_groups'] = list(enzyme_data[r_string].keys())

    activity_data = list(enzyme_data[r_string].values())
    r_dict['value_inner_circle'] = [val[0] for val in activity_data]
    r_dict['value_outer_circle'] = [val[1] for val in activity_data]

    r_dict['attach_atom_id'] = r_group_map[str(r_num)]

    return r_dict

def get_enzyme_dict(enzyme_name, df, scaffold, activity_col='conversion'):
    """ get the r_group data for a single enzyme """

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
    # r_data will look like - {'R1': {'C': 85, 'Br': 45, ...},
    #                          'R2": ...}
    r_data = get_smiles_and_activity(df, r_groups, activity_col=activity_col)

    return r_data

def make_auto_radial_scope(enzyme_1, enzyme_2, df, scaffold_smi,
                           activity_col='conversion', angle_per_wedge=15):

    # make scaffold_mol
    scaffold = Chem.MolFromSmiles(scaffold_smi)

    # get the r_group data for the two enzymes, save to dict
    enz_1_data = get_enzyme_dict(enzyme_1, df, scaffold, activity_col=activity_col)
    enz_2_data = get_enzyme_dict(enzyme_2, df, scaffold, activity_col=activity_col)
    combined_data = combine_enzyme_dict(enz_1_data, enz_2_data)

    # create radial scope settings
    settings = default_settings
    settings['SMILESSTRING'] = scaffold_smi

    # create radial_scope r_dicts, save as list
    r_dicts = []
    r_group_map = create_r_group_map(scaffold)
    for r_num in r_group_map.keys():
        r_setup_dict = get_r_setup_dict(r_num, combined_data, enzyme_1, enzyme_2, scaffold, r_group_map, angle_per_wedge=angle_per_wedge)
        if r_setup_dict is not None:
            r_dicts.append(r_setup_dict)

    # make radial_scope
    scope_plot = rs(settings, *r_dicts)
    return scope_plot

if __name__ == "__main__":
    # Load the data and filter as in Mengyu's notebook
    df = pd.read_excel('AL data.xlsx')
    df = df.replace("Œ±-amino amination", "alpha-amino amination")
    df = df.replace("Œ≤-amino deamination", "beta-amino amination")
    df = df[df['enzyme_type'] == "PAL"]
    df = df.dropna(subset=['product_1_smiles'])
    PandasTools.AddMoleculeColumnToFrame(df, 'product_1_smiles', 'Product Mol')

    scaffold_smi = 'N[C@@H](CC1=C([*:1])C([*:2])=C([*:3])C=C1)C(O)=O'
    make_auto_radial_scope("EncP_E293M", "EncP_E293M", df, scaffold_smi, activity_col='conversion')

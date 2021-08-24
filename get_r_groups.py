import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem import rdRGroupDecomposition
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.warning')

def r_groups_from_scaffold(scaffold):
    """
    Get the r groups used in a scaffold molecule

    For example..
        - r_groups = [1, 2, 3]
    """

    r_groups = []
    for atom in scaffold.GetAtoms():
        atom_props = atom.GetPropsAsDict()

        if 'molAtomMapNumber' in atom_props:
            r_num = atom_props['molAtomMapNumber']
            r_groups.append(r_num)

    return r_groups

def r_group_decomposition(df, scaffold):
    """ Do r group decomposition and return resulting dataframe """

    # Get mols from product smiles
    mols = [Chem.MolFromSmiles(smi) for smi in pd.unique(df['product_1_smiles'])]

    # Use rdkit RGroupDecompose
    groups, unmatched = rdRGroupDecomposition.RGroupDecompose([scaffold], mols, asSmiles=False, asRows=False)

    # Get a list of the unmatched molecules (they don't match the scaffold)
    unmatched_mols = []
    for i in unmatched:
        unmatched_mols.append(mols[i])

    # Get a list of only the molecules that match the scaffold
    mols_match_scaffold = []
    for i, mol in enumerate(mols):
        if i not in unmatched:
            mols_match_scaffold.append(mol)

    # Use RGroupDecompositionToFrame to get dataframe of R group decomposition
    r_group_df = PandasTools.RGroupDecompositionToFrame(groups, mols_match_scaffold, include_core=False,
                                                        redraw_sidechains=False)

    return r_group_df

def make_H_r_groups_none(r_group_df, r_groups):
    """
    To make the r_group_df easier to interpret,
    make any R_group which is just H into None.

    r_group comes from r_groups_from_scaffold
    For example - r_group = [1, 2, ...]

    """

    def remove_H(row, r_groups):

        for num in r_groups:
            r_string = f'R{num}'
            mol = row[r_string]
            smi = Chem.MolToSmiles(mol)

            if smi == f'[H][*:{num}]':
                row[r_string] = None

        return row

    r_group_df.apply(remove_H, args=[r_groups], axis=1)

    return r_group_df

def r_group_df_to_smiles(r_group_df, r_groups):
    """
    Convert the mols in r_group_df to smiles strings
    """

    def show_smiles(mol):
        if mol == None:
            return None
        return Chem.MolToSmiles(mol)

    r_group_df["product_1_smiles"] = r_group_df["Mol"].apply(show_smiles)

    for num in r_groups:
        r_group_df[f"R{num}_smiles"] = r_group_df[f"R{num}"].apply(show_smiles)

    return r_group_df

def split_single_and_multi_r_groups(r_group_df, r_groups):
    """
    Split r_group_df into two dataframes
    - single_r_df - where only 1 r_group is not none
    - multi_r_df - where multiple r_groups are not none
    """

    def has_multiple_rs(row, r_groups):
        '''takes a row and returns True if multiple R groups (and not all H)'''

        count_r = 0
        for num in r_groups:
            r_string = f"R{num}"

            try:
                smi = Chem.MolToSmiles(row[r_string])
            except:
                smi = None

            # if r is not None, or is hydrogen, then add 1 to count of r groups
            if smi != None and smi != f'[H][*:{num}]':
                count_r += 1

        if count_r > 1:
            return True
        else:
            return False

    multiple_r = []
    for index, row in r_group_df.iterrows():
        r = has_multiple_rs(row, r_groups)
        multiple_r.append(r)

    r_group_df['multiple_r'] = multiple_r
    multi_r_df = r_group_df[r_group_df['multiple_r'] == True]
    single_r_df = r_group_df[r_group_df['multiple_r'] == False]

    return single_r_df, multi_r_df

def get_r_df(df, scaffold):
    r_groups = r_groups_from_scaffold(scaffold)
    r_group_df = r_group_decomposition(df, scaffold)
    r_group_df = make_H_r_groups_none(r_group_df, r_groups)
    r_group_df = r_group_df_to_smiles(r_group_df, r_groups)
    single_df, multi_df = split_single_and_multi_r_groups(r_group_df, r_groups)

    return single_df




if __name__ == "__main__":
    # Use __name__ == "__main__' to test a module
    # This code won't run if this is imported, but it will run if this exact .py file run.
    # Great for testing a module..

    # Load the data and filter as in Mengyu's notebook
    df = pd.read_excel('AL data.xlsx')
    df = df.replace("Œ±-amino amination", "alpha-amino amination")
    df = df.replace("Œ≤-amino deamination", "beta-amino amination")
    df = df[df['enzyme_type'] == "PAL"]
    df = df.dropna(subset=['product_1_smiles'])
    PandasTools.AddMoleculeColumnToFrame(df, 'product_1_smiles', 'Product Mol')

    df = df[df['enzyme_name'] == "EncP_E293M"]
    df = df.dropna(subset=['conversion'])
    df = df.drop_duplicates(subset=['product_1_smiles'], keep=False)

    # Test the r group functions
    scaffold = Chem.MolFromSmiles('N[C@@H](CC1=C([*:1])C([*:2])=C([*:3])C([*:4])=C1[*:5])C(O)=O')
    r_groups = r_groups_from_scaffold(scaffold)
    r_group_df = r_group_decomposition(df, scaffold)
    r_group_df = make_H_r_groups_none(r_group_df, r_groups)
    r_group_df = r_group_df_to_smiles(r_group_df, r_groups)
    single_df, multi_df = split_single_and_multi_r_groups(r_group_df, r_groups)

    single_df = get_r_df(df, scaffold)

    print(single_df)

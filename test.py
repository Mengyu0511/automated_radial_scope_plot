import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.warning')
from get_r_groups import get_r_df

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

    single_df = get_r_df(df, scaffold)

    print(single_df)

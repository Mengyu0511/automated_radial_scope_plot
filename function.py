from radialscope import RadialScope as rs
from radialscope import draw_with_indeces
from IPython.display import SVG
from IPython import display
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from rdkit.Chem import PandasTools
IPythonConsole.ipython_useSVG=True
from rdkit.Chem import rdRGroupDecomposition
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.warning')
import math
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import Geometry
import rdkit

def get_repeat_df(dataframe):
    """Ensures the dataframe only has one entry per product, by removing any duplicates"""

    drop_duplicates_df = dataframe.drop_duplicates(subset=['product_1_smiles'], keep=False)

    return drop_duplicates_df

def get_multi_r_df(scaffold, dataframe):
    def groups_to_df(groups,mols,include_core=False,redraw_sidechains=False):
        """ add the molecule to the dataframe 
        """
        cols = ['Mol']+list(groups.keys())
        if redraw_sidechains:
            for k,vl in groups.items():
                if k=='Core':
                    continue
                for i,v in enumerate(vl):
                    vl[i] = Chem.RemoveHs(v)
                    rdDepictor.Compute2DCoords(vl[i])

    
        if not include_core:
            cols.remove('Core')
            del groups['Core']
        groups['Mol'] = mols
        frame = pd.DataFrame(groups,columns=cols)
        PandasTools.RenderImagesInAllDataFrames(images=True)
        return frame
    # creat r_group_df
    mols = [Chem.MolFromSmiles(smi) for smi in pd.unique(dataframe['product_1_smiles'])]
    Draw.MolsToGridImage(mols,molsPerRow=4)
    groups, unmatched = rdRGroupDecomposition.RGroupDecompose([scaffold],mols,asSmiles=False,asRows=False)
    unmatched_mols = []
    for i in unmatched:
        unmatched_mols.append(mols[i])
    
    Draw.MolsToGridImage(unmatched_mols,molsPerRow=4)
    mols_minus_unmatched = []
    for i, mol in enumerate(mols):
        if i not in unmatched:
            mols_minus_unmatched.append(mol)
    r_group_df = PandasTools.RGroupDecompositionToFrame(groups, mols_minus_unmatched, include_core=False, redraw_sidechains=False)
    r_group_df
    def remove_H(mol, r_number=1):
        smi = Chem.MolToSmiles(mol)
        if smi == f'[H][*:{r_number}]':
            return None
        else:
            return mol
    r_group_df['R1'] = r_group_df['R1'].apply(remove_H, r_number=1)
    r_group_df['R2'] = r_group_df['R2'].apply(remove_H, r_number=2)
    r_group_df['R3'] = r_group_df['R3'].apply(remove_H, r_number=3)
    r_group_df['R4'] = r_group_df['R4'].apply(remove_H, r_number=4)
    r_group_df['R5'] = r_group_df['R5'].apply(remove_H, r_number=5)
    r_group_df
    def show_smiles(mol):
        if mol== None:
            return None
        return Chem.MolToSmiles(mol)
    r_group_df["product_1_smiles"] = r_group_df["Mol"].apply(show_smiles)
    r_group_df["R1_smiles"] = r_group_df["R1"].apply(show_smiles)
    r_group_df["R2_smiles"] = r_group_df["R2"].apply(show_smiles)
    r_group_df["R3_smiles"] = r_group_df["R3"].apply(show_smiles)
    r_group_df["R4_smiles"] = r_group_df["R4"].apply(show_smiles)
    r_group_df["R5_smiles"] = r_group_df["R5"].apply(show_smiles)
    r_group_df
    new_df = pd.merge(dataframe,r_group_df,on='product_1_smiles')
    # get multi_r_df and single_r_df

    def has_multiple_rs(row, column_names):
        '''takes a row and returns True if multiple R groups'''
    
        count_r = 0
        for col in column_names:
            if row[col] != None and row[col] != '': 
                count_r += 1
    
        if count_r > 1:
            return True
        else:
            return False
    multiple_r = []
    for index, row in new_df.iterrows():
        r = has_multiple_rs(row, ['R1', 'R2', 'R3', 'R4', 'R5'])
        multiple_r.append(r)

    new_df['multiple_r'] = multiple_r
    multi_r_df = new_df[new_df['multiple_r']==True]
    return multi_r_df

def get_R_list(scaffold,dataframe,column,r_number):
    # 
    def groups_to_df(groups,mols,include_core=False,redraw_sidechains=False):
        """ add the molecule to the dataframe 
        """
        cols = ['Mol']+list(groups.keys())
        if redraw_sidechains:
            for k,vl in groups.items():
                if k=='Core':
                    continue
                for i,v in enumerate(vl):
                    vl[i] = Chem.RemoveHs(v)
                    rdDepictor.Compute2DCoords(vl[i])

    
        if not include_core:
            cols.remove('Core')
            del groups['Core']
        groups['Mol'] = mols
        frame = pd.DataFrame(groups,columns=cols)
        PandasTools.RenderImagesInAllDataFrames(images=True)
        return frame
    # creat r_group_df

    # ---- This should be its own function
    mols = [Chem.MolFromSmiles(smi) for smi in pd.unique(dataframe['product_1_smiles'])]
    Draw.MolsToGridImage(mols,molsPerRow=4)
    groups, unmatched = rdRGroupDecomposition.RGroupDecompose([scaffold],mols,asSmiles=False,asRows=False)
    unmatched_mols = []
    for i in unmatched:
        unmatched_mols.append(mols[i])
    
    Draw.MolsToGridImage(unmatched_mols,molsPerRow=4)
    mols_minus_unmatched = []
    for i, mol in enumerate(mols):
        if i not in unmatched:
            mols_minus_unmatched.append(mol)
    r_group_df = PandasTools.RGroupDecompositionToFrame(groups, mols_minus_unmatched, include_core=False, redraw_sidechains=False)
    r_group_df

    # ---- To Here

    # ---- This should be a function, but make it work for any number of R groups
    def remove_H(mol, r_number=1):
        smi = Chem.MolToSmiles(mol)
        if smi == f'[H][*:{r_number}]':
            return None
        else:
            return mol
    r_group_df['R1'] = r_group_df['R1'].apply(remove_H, r_number=1)
    r_group_df['R2'] = r_group_df['R2'].apply(remove_H, r_number=2)
    r_group_df['R3'] = r_group_df['R3'].apply(remove_H, r_number=3)
    r_group_df['R4'] = r_group_df['R4'].apply(remove_H, r_number=4)
    r_group_df['R5'] = r_group_df['R5'].apply(remove_H, r_number=5)
    r_group_df

    # ---- To here


    def show_smiles(mol):
        if mol== None:
            return None
        return Chem.MolToSmiles(mol)
    r_group_df["product_1_smiles"] = r_group_df["Mol"].apply(show_smiles)
    r_group_df["R1_smiles"] = r_group_df["R1"].apply(show_smiles)
    r_group_df["R2_smiles"] = r_group_df["R2"].apply(show_smiles)
    r_group_df["R3_smiles"] = r_group_df["R3"].apply(show_smiles)
    r_group_df["R4_smiles"] = r_group_df["R4"].apply(show_smiles)
    r_group_df["R5_smiles"] = r_group_df["R5"].apply(show_smiles)
    r_group_df
    new_df = pd.merge(dataframe,r_group_df,on='product_1_smiles')
    # get multi_r_df and single_r_df


    # ---- This can be its own function
    def has_multiple_rs(row, column_names):
        '''takes a row and returns True if multiple R groups'''
    
        count_r = 0
        for col in column_names:
            if row[col] != None and row[col] != '': 
                count_r += 1
    
        if count_r > 1:
            return True
        else:
            return False
    multiple_r = []
    for index, row in new_df.iterrows():
        r = has_multiple_rs(row, ['R1', 'R2', 'R3', 'R4', 'R5'])
        multiple_r.append(r)

    new_df['multiple_r'] = multiple_r
    multi_r_df = new_df[new_df['multiple_r']==True]
    single_r_df = new_df[new_df['multiple_r']==False]

    # ---- To here


    # ---- This code is now to getting the smiles for a particular R group
    # ---- Can be its own function

    # remove the repeat product
    single_R_df = single_r_df.drop_duplicates("product_1_smiles")

    single_r_group_df=pd.DataFrame(single_R_df,columns=["Mol","R1","R2","R3","R4","R5"])
    
    R_group_df=pd.DataFrame(single_r_group_df,columns=["Mol",column]).dropna()

    # Haven't we already done this above?
    # Change them to the smiles strings
    R_group_df[column+"_smiles"] = R_group_df[column].apply(show_smiles)

    R_list=R_group_df[column+"_smiles"].tolist()
    for i, sml in enumerate(R_list):
        if sml.endswith(f"[*:{r_number}]"):
             R_list[i] = R_list[i][:-5]

    return  R_list

def get_conversion(scaffold,dataframe,column,r_number):
    # 
    def groups_to_df(groups,mols,include_core=False,redraw_sidechains=False):
        """ add the molecule to the dataframe 
        """
        cols = ['Mol']+list(groups.keys())
        if redraw_sidechains:
            for k,vl in groups.items():
                if k=='Core':
                    continue
                for i,v in enumerate(vl):
                    vl[i] = Chem.RemoveHs(v)
                    rdDepictor.Compute2DCoords(vl[i])

    
        if not include_core:
            cols.remove('Core')
            del groups['Core']
        groups['Mol'] = mols
        frame = pd.DataFrame(groups,columns=cols)
        PandasTools.RenderImagesInAllDataFrames(images=True)
        return frame
    # creat r_group_df
    mols = [Chem.MolFromSmiles(smi) for smi in pd.unique(dataframe['product_1_smiles'])]
    Draw.MolsToGridImage(mols,molsPerRow=4)
    groups, unmatched = rdRGroupDecomposition.RGroupDecompose([scaffold],mols,asSmiles=False,asRows=False)
    unmatched_mols = []
    for i in unmatched:
        unmatched_mols.append(mols[i])
    
    Draw.MolsToGridImage(unmatched_mols,molsPerRow=4)
    mols_minus_unmatched = []
    for i, mol in enumerate(mols):
        if i not in unmatched:
            mols_minus_unmatched.append(mol)
    r_group_df = PandasTools.RGroupDecompositionToFrame(groups, mols_minus_unmatched, include_core=False, redraw_sidechains=False)
    r_group_df
    def remove_H(mol, r_number=1):
        smi = Chem.MolToSmiles(mol)
        if smi == f'[H][*:{r_number}]':
            return None
        else:
            return mol
    r_group_df['R1'] = r_group_df['R1'].apply(remove_H, r_number=1)
    r_group_df['R2'] = r_group_df['R2'].apply(remove_H, r_number=2)
    r_group_df['R3'] = r_group_df['R3'].apply(remove_H, r_number=3)
    r_group_df['R4'] = r_group_df['R4'].apply(remove_H, r_number=4)
    r_group_df['R5'] = r_group_df['R5'].apply(remove_H, r_number=5)
    r_group_df
    def show_smiles(mol):
        if mol== None:
            return None
        return Chem.MolToSmiles(mol)
    r_group_df["product_1_smiles"] = r_group_df["Mol"].apply(show_smiles)
    r_group_df["R1_smiles"] = r_group_df["R1"].apply(show_smiles)
    r_group_df["R2_smiles"] = r_group_df["R2"].apply(show_smiles)
    r_group_df["R3_smiles"] = r_group_df["R3"].apply(show_smiles)
    r_group_df["R4_smiles"] = r_group_df["R4"].apply(show_smiles)
    r_group_df["R5_smiles"] = r_group_df["R5"].apply(show_smiles)
    r_group_df
    new_df = pd.merge(dataframe,r_group_df,on='product_1_smiles')
    # get multi_r_df and single_r_df

    def has_multiple_rs(row, column_names):
        '''takes a row and returns True if multiple R groups'''
    
        count_r = 0
        for col in column_names:
            if row[col] != None and row[col] != '': 
                count_r += 1
    
        if count_r > 1:
            return True
        else:
            return False
    multiple_r = []
    for index, row in new_df.iterrows():
        r = has_multiple_rs(row, ['R1', 'R2', 'R3', 'R4', 'R5'])
        multiple_r.append(r)

    new_df['multiple_r'] = multiple_r
    multi_r_df = new_df[new_df['multiple_r']==True]
    single_r_df = new_df[new_df['multiple_r']==False]
    # remove the repeat product
    single_R_df = single_r_df.drop_duplicates("product_1_smiles")
    single_R_group_df=pd.DataFrame(single_R_df,columns=["Mol","R1","R2","R3","R4","R5"])
    
    R_group_df=pd.DataFrame(single_R_group_df,columns=["Mol",column]).dropna()
    R_group_df["product_1_smiles"] = R_group_df["Mol"].apply(show_smiles)


    # We only want this last part - everything else is repeated.

    conversion_df = pd.merge(single_R_df,R_group_df,on='product_1_smiles')
    conversion_list=conversion_df["conversion"].tolist()
    new_conversions = [int(c) for c in conversion_list]

    return  new_conversions

def get_bond_angle(mol, r_group, debug=False):
    def show_atom_ids(mol, label):
        for atom in mol.GetAtoms():
                atom.SetProp(label, str(atom.GetIdx()))
        return mol

    def create_r_group_map(mol):
        """
        R group ids are not the same as atom ids
        Need a map to go from the R group number to the atom number
        This function creates that map as a dictionary.
        (The keys of dictionaries can't be numbers, so convert to str first for key)
        """
        r_group_map = {}
    
        for atom in mol.GetAtoms():
            atom_props = atom.GetPropsAsDict()
            atom_id = atom.GetIdx()
            
            # if atom has prop 'molAtomMapNumber', then this is the r group label
            if 'molAtomMapNumber' in atom_props:
                r_group_map[str(atom_props['molAtomMapNumber'])] = atom_id
            
        return r_group_map

    
    
    """
    This function will get the r group atom, find its neighbour, 
    and calculate the bond angle between them in degrees.
    
    (I haven't tested whether the angle is actually correct yet..)
    (Also, I'm not 100% what angle 0 degrees corresponds to, and if this matches the radial scope plot)
    """
    
    # create necessary variables
    r_group_map = create_r_group_map(mol)
    dm = Draw.PrepareMolForDrawing(mol)
    
    # get r group atom id
    atom_id = r_group_map[str(r_group)]
    atom = mol.GetAtomWithIdx(atom_id)
        
    # find the coordinates of the atom
    atom_pos = Geometry.Point2D(dm.GetConformer().GetAtomPosition(atom_id))
    
    # get the atoms neighbours, should only be one
    neighbours = atom.GetNeighbors()
    if len(neighbours) != 1:
        print(f"Warning - atom {atom_num} has {len(neighbours)} neighbours")

    # get the coordinates of the neighbour
    first_neighbour = neighbours[0]
    first_neighbour_id = first_neighbour.GetIdx()
    neighbour_pos = Geometry.Point2D(dm.GetConformer().GetAtomPosition(first_neighbour_id))
    
    # calculate the vector for the bond between the atoms, using the coordinates
    direction_vector = Geometry.Point2D.DirectionVector(atom_pos, neighbour_pos)
    
    # angle of vector - taken from https://stackoverflow.com/questions/6247153/angle-from-2d-unit-vector
    # not sure if this is correct..  seems ok?
    angle = math.atan2(direction_vector.x, direction_vector.y)*180/math.pi 
    
    # print commands for debugging if its not working
    if debug == True:
        print(f"R group: {r_group}, Atom_ID: {atom_id}")
        print(f"X: {atom_pos.x}, Y: {atom_pos.y}")
        print(f"Number of neighbours = {len(neighbours)}")
        print(f"Neighbour_ID: {first_neighbour_id}, Neighbour x: {neighbour_pos.x}, Neighbour y: {neighbour_pos.y}")
        print(f"Angle = {angle}")
    
    return angle
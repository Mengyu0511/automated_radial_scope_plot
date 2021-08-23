from rdkit.Chem import Draw
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.warning')
import math
from rdkit import Geometry

def show_atom_ids(mol, label):
    """
    Use to display a molecule with its atom ids
    """
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

def get_bond_angle(mol, r_group, debug=False):
    """
    This function will get the r group atom, find its neighbour,
    and calculate the bond angle between them in degrees.
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
        print(f"Warning - atom {atom_id} has {len(neighbours)} neighbours")

    # get the coordinates of the neighbour
    first_neighbour = neighbours[0]
    first_neighbour_id = first_neighbour.GetIdx()
    neighbour_pos = Geometry.Point2D(dm.GetConformer().GetAtomPosition(first_neighbour_id))

    # calculate the vector for the bond between the atoms, using the coordinates
    direction_vector = Geometry.Point2D.DirectionVector(atom_pos, neighbour_pos)

    # angle of vector - taken from https://stackoverflow.com/questions/6247153/angle-from-2d-unit-vector
    angle = math.atan2(direction_vector.x, direction_vector.y ) * 180 / math.pi

    # print commands for debugging if its not working
    if debug == True:
        print(f"R group: {r_group}, Atom_ID: {atom_id}")
        print(f"X: {atom_pos.x}, Y: {atom_pos.y}")
        print(f"Number of neighbours = {len(neighbours)}")
        print(f"Neighbour_ID: {first_neighbour_id}, Neighbour x: {neighbour_pos.x}, Neighbour y: {neighbour_pos.y}")
        print(f"Angle = {angle}")

    return angle

def calculate_r_scope_angle(bond_angle, num_wedges, angle_per_wedge):

    r_scope_angle = 270 - bond_angle + angle_per_wedge * num_wedges / 2
    return r_scope_angle

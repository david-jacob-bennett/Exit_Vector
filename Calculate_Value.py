from sys import argv
from Bio.PDB import PDBList
from Bio.PDB import MMCIFParser

def pull_pdb_file(pdb_id):
    """
    We need to pull the cif file from the PDB database.
    """

    pdbl = PDBList()

    file_path = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='mmCif')

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(pdb_id, file_path)
    CA_list = []
    for model in structure:
        for chain in model:
            for residue in chain:
                residue_id = residue.get_id()[1]
                if 76 <= residue_id <= 87 and "CA" in residue:
                    ca_coord = residue["CA"].get_coord()
                    CA_list.append(tuple([residue_id, ca_coord]))
    return CA_list


def find_coord_averages(pdb_id):
    ca_list = pull_pdb_file(pdb_id)
    bottom_helix = [ca for ca in ca_list if 76 <= ca[0] <= 81]
    top_helix = [ca for ca in ca_list if 82 <= ca[0] <= 87]
    x_avg_bottom = sum(ca[1][0] for ca in bottom_helix) / len(bottom_helix)
    y_avg_bottom = sum(ca[1][1] for ca in bottom_helix) / len(bottom_helix)
    z_avg_bottom = sum(ca[1][2] for ca in bottom_helix) / len(bottom_helix)

    x_avg_top = sum(ca[1][0] for ca in top_helix) / len(top_helix)
    y_avg_top = sum(ca[1][1] for ca in top_helix) / len(top_helix)
    z_avg_top = sum(ca[1][2] for ca in top_helix) / len(top_helix)
    top_coords = (x_avg_top, y_avg_top, z_avg_top)
    bottom_coords = (x_avg_bottom, y_avg_bottom, z_avg_bottom)
    print (bottom_coords, top_coords)
    return bottom_coords, top_coords

def get_AC_dx_dy_dz(bottom_coords, top_coords):
    dx = top_coords[0] - bottom_coords[0]
    dy = top_coords[1] - bottom_coords[1]
    dz = top_coords[2] - bottom_coords[2]
    return dx, dy, dz

def get_AB_dx_dy_dz(bottom_coords, Top_axis):
    dx = bottom_coords[0] - Top_axis[0]
    dy = bottom_coords[1] - Top_axis[1]
    dz = bottom_coords[2] - Top_axis[2]
    return dx, dy, dz

def get_BC_dx_dy_dz(top_coords, top_axis):
    dx = top_coords[0] - top_axis[0]
    dy = top_coords[1] - top_axis[1]
    dz = top_coords[2] - top_axis[2]
    return dx, dy, dz

def output_top_axis(bottom_coords):
    top_axis_x = bottom_coords[0]
    top_axis_y = bottom_coords[1]
    top_axis_z = 10
    return top_axis_x, top_axis_y, top_axis_z

def find_length(AC, AB, BC):
    length_AC = (AC[0]**2 + AC[1]**2 + AC[2]**2)**0.5
    length_AB = (AB[0]**2 + AB[1]**2 + AB[2]**2)**0.5
    length_BC = (BC[0]**2 + BC[1]**2 + BC[2]**2)**0.5
    return length_AC, length_AB, length_BC

def find_exit_vector(side_lengths):
    length_AC, length_AB, length_BC = side_lengths
    # Using the law of cosines to find the angle at point B
    from math import acos, degrees
    angle_ABC = degrees(acos((length_AB**2 + length_BC**2 - length_AC**2) / (2 * length_AB * length_BC)))
    return angle_ABC

def main():
    bottom_coords, top_coords = find_coord_averages(argv[1])
    top_axis = output_top_axis(bottom_coords)
    AC = get_AC_dx_dy_dz(bottom_coords, top_coords)
    AB = get_AB_dx_dy_dz(bottom_coords, top_axis)
    BC = get_BC_dx_dy_dz(top_coords, top_axis)
    side_lengths = find_length(AC, AB, BC)
    angle_ABC = find_exit_vector(side_lengths)
    print(angle_ABC)

if __name__ == "__main__":
    main()
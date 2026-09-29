import pyrosetta
from pyrosetta import * 
from pyrosetta.rosetta.core.scoring.dssp import Dssp
import numpy as np
import sys
from sys import argv
import rosetta
import os
import subprocess

pyrosetta.init(
    "-crystal_refine",
    "-cryst:refinable_lattice",
    "-cryst:interaction_shell" ,
    "-score_symm_complex",
    "-mute all"
    
)
# not currently using. Might in the future
def generate_sym_file(sym_file_name):
    if not os.path.exists(sym_file_name):
        rosetta_path = os.environ.get("ROSETTA3", "/home/benchmark/rosetta")
        perl_script = os.path.join(rosetta_path, "source/src/ apps/public/symmetry/make_symmdef_file.pl")
        cmd = f"perl {perl_script} -m CRYST -p 9dp8.pdb > {sym_file_name}"
        subprocess.run(cmd, shell=True, check=True)
# errors with pdb being appended too many times. Starting to work though
def find_linker_s_e(pose, linker):
   
    rosetta.basic.options.set_real_option("cryst:interaction_shell", 12.0)
    updated_pose = pyrosetta.pose_from_file(pose + ".pdb")

    
    symm_mover = rosetta.protocols.symmetry.SetupForSymmetryMover(updated_pose.pdb_info().name() + ".sym")
    symm_mover.apply(updated_pose)

    
    pose_seq = updated_pose.sequence()
    start_idx = pose_seq.find(linker)
    if start_idx == -1:
        raise ValueError("Linker sequence not found in the pose")
    start_res = start_idx + 1
    end_res = start_res + len(linker) - 1
    return start_res, end_res
# code hasn't gotten this far, don't know if it works yet
def get_exit_vector(pose, start_res, end_res):
    coords = np.array([pose.residue(i).xyz("CA") for i in range(start_res, end_res + 1)])
    eigenvalues, eigenvectors = np.linalg.eigh(np.dot((coords - coords.mean(axis=0)).T, coords - coords.mean(axis=0)))

    exit_vector = eigenvectors[:,-1]
    if np.dot(exit_vector, coords[-1] - coords[0]) < 0:
        exit_vector = -exit_vector
    exit_vector /= np.linalg.norm(exit_vector)
    return exit_vector
    


def main():
    start_res, end_res = find_linker_s_e(argv[1], argv[2])
    get_exit_vector(argv[1], start_res, end_res)

if __name__ == "__main__":
    main()
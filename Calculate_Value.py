import sys
import numpy as np
import pyrosetta

pyrosetta.init("-mute all")

def ca_xyz(pose, i):
    v = pose.residue(i).xyz("CA")
    return np.array([v.x, v.y, v.z])

def find_linker(pose, linker):
    idx = pose.sequence().find(linker)
    if idx == -1:
        raise ValueError("Linker sequence not found in the pose")
    return idx + 1, idx + len(linker)

def exit_vector(pose, start, end):
    coords = np.array([ca_xyz(pose, i) for i in range(start, end + 1)])
    centered = coords - coords.mean(axis=0)
    _, vecs = np.linalg.eigh(centered.T @ centered)
    axis = vecs[:, -1]                      # already unit length
    if np.dot(axis, coords[-1] - coords[0]) < 0:
        axis = -axis                        # point N -> C
    return axis

def outward_angle(pose, start, end, axis, n_base=4):
    all_ca = np.array([ca_xyz(pose, i) for i in range(1, pose.total_residue() + 1)])
    center = all_ca.mean(axis=0)
    base = np.mean([ca_xyz(pose, i) for i in range(start, start + n_base)], axis=0)
    out = base - center
    out /= np.linalg.norm(out)
    return np.degrees(np.arccos(np.clip(np.dot(axis, out), -1.0, 1.0)))

def main():
    if len(sys.argv) < 3:
        sys.exit("Usage: python Calculate_Value.py <pose_prefix> <linker_seq>")
    pose = pyrosetta.pose_from_file(sys.argv[1] + ".pdb")
    start, end = find_linker(pose, sys.argv[2])
    axis = exit_vector(pose, start, end)
    print("Helix exit vector:", axis)
    print(f"Angle to outward vector: {outward_angle(pose, start, end, axis):.2f}°")

if __name__ == "__main__":
    main()

# import sys
# import numpy as np
# import pyrosetta

# pyrosetta.init("-mute all")

# HELIX_LEN = 0   # extra residues before the linker to include in the axis fit (0 = linker only)


# def ca_xyz(pose, i):
#     v = pose.residue(i).xyz("CA")
#     return np.array([v.x, v.y, v.z])


# def find_linker(pose, linker):
#     seq = pose.sequence()
#     n = seq.count(linker)
#     if n == 0:
#         raise ValueError(f"Linker {linker} not found in pose sequence")
#     if n > 1:
#         print(f"WARNING: {linker} appears {n} times; using the first match")
#     idx = seq.find(linker)
#     start, end = idx + 1, idx + len(linker)
#     print(f"Linker {linker}: residues {start}-{end}")
#     return start, end


# def exit_vector(pose, start, end):
#     coords = np.array([ca_xyz(pose, i) for i in range(start, end + 1)])
#     centered = coords - coords.mean(axis=0)
#     _, vecs = np.linalg.eigh(centered.T @ centered)
#     axis = vecs[:, -1]
#     if np.dot(axis, coords[-1] - coords[0]) < 0:
#         axis = -axis                         # point N -> C
#     return axis


# def outward_angle(pose, start, end, axis, n_base=4):
#     all_ca = np.array([ca_xyz(pose, i) for i in range(1, pose.total_residue() + 1)])
#     center = all_ca.mean(axis=0)
#     n = min(n_base, end - start + 1)
#     base = np.mean([ca_xyz(pose, i) for i in range(start, start + n)], axis=0)
#     out = base - center
#     out /= np.linalg.norm(out)
#     return np.degrees(np.arccos(np.clip(np.dot(axis, out), -1.0, 1.0)))


# def main():
#     if len(sys.argv) < 3:
#         sys.exit("Usage: python Calculate_Value.py <pose_prefix> <linker_seq>")
#     pose = pyrosetta.pose_from_file(sys.argv[1] + ".pdb")
#     link_start, link_end = find_linker(pose, sys.argv[2])

#     win_start = max(1, link_start - HELIX_LEN)
#     axis = exit_vector(pose, win_start, link_end)
#     angle = outward_angle(pose, win_start, link_end, axis)

#     print(f"Axis window: residues {win_start}-{link_end}")
#     print("Helix exit vector:", axis)
#     print(f"Angle to outward vector: {angle:.2f}°")


# if __name__ == "__main__":
#     main()
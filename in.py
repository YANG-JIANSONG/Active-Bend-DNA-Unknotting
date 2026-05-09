from lammps import lammps
import numpy as np
import random
import ctypes
import os
import math
import pythonknot.alexander_poly as alexander_poly

# Function to save XYZ trajectory
def save_xyz_traj(file, X: np.ndarray):
    """Appends atom positions to an XYZ trajectory file."""
    with open(file, 'a') as f:
        for x in X:
            f.write(str(x.shape[0]) + '\n')
            f.write('\n')
            for p in x:
                f.write('1\t{:f} {:f} {:f}'.format(p[0], p[1], p[2]) + '\n')

def plane_from_points(points):
    """Fits a plane using three points; returns the normal vector and a point on the plane."""
    p0 = points[0]
    v1 = points[1] - points[0]
    v2 = points[2] - points[0]
    normal = np.cross(v1, v2)
    normal = normal / np.linalg.norm(normal)
    return normal, p0

def line_plane_intersection(p1, p2, plane_normal, plane_point):
    """Calculates the intersection of a line segment and a plane."""
    u = p2 - p1
    denom = np.dot(plane_normal, u)
    if abs(denom) < 1e-12:  # Parallel case
        return None
    t = np.dot(plane_normal, plane_point - p1) / denom
    if 0 <= t <= 1:  # Intersection point is within the segment
        return p1 + t * u
    return None

def project_to_plane(points):
    """Fits a 3D point set to a plane and returns coordinates projected onto a 2D plane."""
    points = np.array(points)
    centroid = points.mean(axis=0)
    u, s, vh = np.linalg.svd(points - centroid)
    normal = vh[2, :]  # Normal vector of the plane
    
    # Build local coordinate system for the plane
    x_axis = vh[0, :]
    y_axis = vh[1, :]
    
    # Project to local coordinate system
    proj_2d = np.dot(points - centroid, np.vstack((x_axis, y_axis)).T)
    return proj_2d, centroid, x_axis, y_axis

def is_point_in_polygon(x, y, poly):
    """Ray-casting algorithm to determine if a point is inside a 2D polygon."""
    n = len(poly)
    inside = False
    for i in range(n):
        j = (i + 1) % n
        xi, yi = poly[i]
        xj, yj = poly[j]
        if ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / (yj - yi + 1e-15) + xi):
            inside = not inside
    return inside

def compute_dis_and_theta(data, i, j):
    """Computes distance and angle between segments for topological analysis."""
    p1, p2, p3 = data[i-1], data[i], data[i+1]
    plane_point = p2
    normal = get_plane_normal(p1, p2, p3)

    q1, q2 = data[j-1], data[j]
    direction = q2 - q1

    intersection = line_plane_intersection(q1, direction, plane_point, normal)
    if intersection is None:
        return None, None

    dis = np.linalg.norm(intersection - p2)

    cos_theta = np.dot(direction, normal) / (np.linalg.norm(direction) * np.linalg.norm(normal))
    angle_between = np.arccos(np.clip(np.abs(cos_theta), -1.0, 1.0))
    theta_deg = 90.0 - np.degrees(angle_between)

    # Check if intersection is within the sector defined by p1-p2-p3
    v1 = p1 - p2
    v2 = p3 - p2
    if not is_within_sector(p2, v1, v2, intersection):
        dis = -dis  # Intersection point is outside the sector

    return dis, theta_deg

def is_within_sector(center, v1, v2, target):
    """Checks if a target point lies within the angular sector formed by vectors v1 and v2."""
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)
    vt = (target - center)
    vt = vt / np.linalg.norm(vt)

    cross_ref = np.cross(v1, v2)
    if np.dot(np.cross(v1, vt), cross_ref) < 0:
        return False
    if np.dot(np.cross(v2, vt), cross_ref) > 0:
        return False
    return True

def distance(x1: np.ndarray, x2):
    return np.sqrt(np.sum((x1-x2)**2))

def calculate_plane_angle(p1, p2, p3, q1, q2, q3):
    """Calculates the angle in degrees between two planes defined by three points each."""
    def calculate_normal(a, b, c):
        v1 = np.array(b) - np.array(a)
        v2 = np.array(c) - np.array(a)
        normal = np.cross(v1, v2)
        return normal / np.linalg.norm(normal)
    
    normal1 = calculate_normal(p1, p2, p3)
    normal2 = calculate_normal(q1, q2, q3)
    
    dot_product = np.dot(normal1, normal2)
    dot_product = np.clip(dot_product, -1.0, 1.0)
    angle_radians = np.arccos(dot_product)
    return np.degrees(angle_radians)

def calculate_curvature(points: np.ndarray):
    """Calculates the bending angle for each monomer along the chain."""
    curvature = 0.0
    for index, x in enumerate(points):
        if index > 0 and index < (points.shape[0]-1):
            a = points[index-1] - points[index]
            b = points[index+1] - points[index]
            bend_angle = np.pi - vectorsangle(a, b)
            curvature = bend_angle
    return curvature

def get_plane_normal(p1, p2, p3):
    v1 = p2 - p1
    v2 = p3 - p1
    normal = np.cross(v1, v2)
    return normal / np.linalg.norm(normal)

def bendenergy_calc(knotcore, lp):
    """Calculates bending energy based on persistence length and curvature."""
    curvature = calculate_curvature(knotcore)
    E_bend = 0.5 * lp * (curvature**2)
    return E_bend

def vectorsangle(a, b):
    dot_product = np.dot(a, b)
    costheta = dot_product / (np.linalg.norm(a) * np.linalg.norm(b))
    theta = math.acos(np.clip(costheta, -1.0, 1.0))
    return theta

def is_point_in_angle(P, data):
    """Checks if point P is within the 2D angle formed by points A-B-C."""
    A, B, C = data
    pts = [A, B, C]
    O = np.array(pts[1])  # Vertex of the angle
    P = np.array(P)
    
    j, k = [x for x in range(3) if x != 1]
    v1 = np.array(pts[j]) - O
    v2 = np.array(pts[k]) - O
    vp = P - O
    
    cross1 = np.cross(v1, vp)
    cross2 = np.cross(vp, v2)
    cross_angle = np.cross(v1, v2)
    
    # Check sign of cross products to determine if vp is between v1 and v2
    inside = (cross1 * cross_angle >= 0) and (cross2 * cross_angle >= 0)
    return inside

def hugging(bend, T):
    """Checks the 'hugging' condition for potential strand passage."""
    if (T[0][2] < -1 and T[2][2] > 1) or (T[0][2] > 1 and T[2][2] < -1):
        bend = np.array(bend)
        bend2d = bend[:, :2]
        T = np.array(T)
        T2d = T[:, :2]
        if is_point_in_angle(T2d[1], [bend2d[0], bend2d[1], bend2d[2]]):
            return True
    return False

def TOPO(lmp):
    """Topological check to find potential crossing sites."""
    natoms = lmp.get_natoms()
    L = natoms
    x = lmp.gather_atoms("x", 1, 3)
    positions = np.reshape(x, (natoms, 3))
    data = positions

    core = 10
    c = 0
    for i in range(2, 3):
        for j in range(i, L - 1):
            if abs(i - j) >= core and abs(i - j) <= L - core:
                dis = distance(data[i], data[j])
                if dis <= 4:
                    if hugging([data[i-2], data[i], data[i+2]], [data[j-2], data[j], data[(j+2)%L]]):
                        c += 1
                        return i, j
    if c == 0:
        return 0, 0

def swap_atoms_old(lmp, index):
    """Attempts to perform a strand passage by swapping atom coordinates if conditions are met."""
    natoms = lmp.get_natoms()
    L = natoms
    x = lmp.gather_atoms("x", 1, 3)
    positions = np.reshape(x, (natoms, 3))
    
    data = positions.copy()
    ori_data = data.copy()
    
    core = 10
    c = 0
    for i in range(2, 3):
        for j in range(i, L - 1):
            if abs(i - j) >= core and abs(i - j) <= L - core:
                dis = distance(data[i], data[j])
                if dis <= 3:
                    ang = calculate_plane_angle(data[(i-1)%L], data[i], data[(i+1)%L], data[(j-1)%L], data[j], data[(j+1)%L])
                    
                    # Perform coordinate swap for the local segments
                    data[[i, j]] = data[[j, i]]
                    data[[i-2, j-2]] = data[[j-2, i-2]]
                    data[[i-1, j-1]] = data[[j-1, i-1]]
                    data[[(i+1)%L, (j+1)%L]] = data[[(j+1)%L, (i+1)%L]]
                    data[[(i+2)%L, (j+2)%L]] = data[[(j+2)%L, (i+2)%L]]
                    
                    knot_core_new = np.reshape(data, ((1, L, 3)))
                    knot_type1 = alexander_poly.calculate_knot_type(knot_core_new, "ring")
                    
                    knot_core_old = np.reshape(ori_data, ((1, L, 3)))
                    knot_type0 = alexander_poly.calculate_knot_type(knot_core_old, "ring")

                    # Log successful unknotting or topological transition events
                    if knot_type0 == ['1'] and knot_type1 != ['1']:
                        print("Transition from Unknot detected")
                        with open('output_index_0921.txt', 'a') as f:
                            print(index, file=f)
                        with open('output_all_save_0921.txt', 'a') as f:
                            print(knot_type0, knot_type1, dis, ang, [i, j], file=f)
                        save_xyz_traj('trj_reco_0921.xyz', knot_core_old)
                        
                    if knot_type0 == ['3_1'] and knot_type1 == ['1']:
                        print("Unknotting (3_1 to 1) detected")
                        with open('output_index_0921.txt', 'a') as f:
                            print(index, file=f)
                        with open('output_all_save_0921.txt', 'a') as f:
                            print(knot_type0, knot_type1, dis, ang, [i, j], file=f)
                        save_xyz_traj('trj_reco_0921.xyz', knot_core_old)
                    
                    c += 1
                    return i, j
    if c == 0:
        return 0, 0

def Run_lammps(lp, seed, L, type, degree):
    """Main simulation loop using LAMMPS Python interface."""
    lmp = lammps(cmdargs=["-screen", "none", "-log", "none"])
    halflp = lp / 2
    
    # Initialize LAMMPS system
    commands = [
        "log none",
        "units lj",
        "boundary f f f",
        "region box block -100 100 -100 100 -50 50",
        "atom_style angle",
        f"read_data lammps_{type}_L{L}_close.data",
        # Group definitions for the active bend site
        "group rigid_atoms id 1 2 3 4 5",
        "velocity rigid_atoms set 0.0 0.0 0.0",
        "fix hold rigid_atoms setforce 0.0 0.0 0.0",
        
        "group not_rigid subtract all rigid_atoms",
        "group atom1 molecule 1",
        "bond_style fene",
        "bond_coeff 1 30. 3. 1.0 1.0",
        "angle_style harmonic",
        f"angle_coeff 1 {halflp} 180.0",
        "special_bonds lj 0.0 1.0 1.0",
        "pair_style lj/cut 1.122462",
        "pair_coeff 1 1 1.0 1.0 1.122462",
        "neighbor 2.0 bin",
        "fix 1 not_rigid nve/limit 0.1",
        f"fix 2 not_rigid langevin 1.0 1.0 2.0 {seed}",
        "fix 3 not_rigid wall/reflect xlo -100 xhi 100 ylo -100 yhi 100 zlo -50 zhi 50",
        "timestep 0.01",
        "thermo_style custom step time temp etotal eangle ke cpu",
        "run 0",
        "thermo 500000",
    ]

    for cmd in commands:
        lmp.command(cmd)

    # Simulation loop with topological sampling
    for step in range(1000000):
        swap_atoms_old(lmp, step)  # Check and perform strand passage
        lmp.command("run 100")     # Propagate dynamics for 100 steps

# Example Run
Run_lammps(lp=17, seed=1, L=300, type=31, degree=0)
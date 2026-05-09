import numpy as np
import rmsd

def bezier_point(P0, P1, P2, t):
    """Calculates a point on a quadratic Bezier curve at parameter t."""
    return (1 - t)**2 * P0 + 2 * (1 - t) * t * P1 + t**2 * P2

def bezier_length(P0, P1, P2, samples=1000):
    """Estimates the total arc length of the Bezier curve via linear sampling."""
    ts = np.linspace(0, 1, samples)
    points = np.array([bezier_point(P0, P1, P2, t) for t in ts])
    diffs = points[1:] - points[:-1]
    seg_lengths = np.linalg.norm(diffs, axis=1)
    return np.sum(seg_lengths), points, ts

def sample_points_by_arc_length(P0, P1, P2, num_points):
    """Samples points along the curve at equal arc-length intervals."""
    length, curve_points, ts = bezier_length(P0, P1, P2)
    seg_lengths = np.linalg.norm(curve_points[1:] - curve_points[:-1], axis=1)
    cum_length = np.concatenate(([0], np.cumsum(seg_lengths)))

    target_lengths = np.linspace(0, length, num_points)
    points = []
    for s in target_lengths:
        idx = np.searchsorted(cum_length, s)
        if idx == 0:
            points.append(curve_points[0])
        elif idx >= len(cum_length):
            points.append(curve_points[-1])
        else:
            t0, t1 = ts[idx - 1], ts[idx]
            l0, l1 = cum_length[idx - 1], cum_length[idx]
            t = t0 + (s - l0) / (l1 - l0) * (t1 - t0)
            pt = bezier_point(P0, P1, P2, t)
            points.append(pt)
    return np.array(points)

def insert_points_triangle_method(p1, p2, num_points, spacing=1.0):
    """Generates an isosceles triangle path between two points to fix segment lengths."""
    p1 = np.array(p1, dtype=float)
    p2 = np.array(p2, dtype=float)
    L = np.linalg.norm(p2 - p1)
    
    total_segments = num_points
    total_length = total_segments * spacing
    if total_length < L:
        raise ValueError("Spacing too small to construct a valid path.")

    # Calculate height of the isosceles triangle
    half_length = total_length / 2
    h = np.sqrt(half_length**2 - (L/2)**2)

    mid = (p1 + p2) / 2
    base_dir = (p2 - p1) / L

    # Define normal vector (defaults to Z-axis unless aligned)
    if np.allclose(base_dir, [0, 0, 1]) or np.allclose(base_dir, [0, 0, -1]):
        normal = np.array([0, 1, 0])
    else:
        normal = np.cross(base_dir, [0, 0, 1])
        normal /= np.linalg.norm(normal)

    apex = mid + h * normal
    n1 = total_segments // 2
    n2 = total_segments - n1

    def interpolate(start, end, n_segments):
        return [start + (end - start) * i / n_segments for i in range(n_segments + 1)]

    left_points = interpolate(p1, apex, n1)
    right_points = interpolate(apex, p2, n2)[1:]

    points = left_points + right_points
    return points[1:-1]

def check_distance(coords, tol=1e-2):
    """Validates if bond lengths between consecutive atoms are approximately 1.0."""
    bad = []
    for i in range(len(coords)):
        j = (i + 1) % len(coords)
        d = np.linalg.norm(coords[i] - coords[j])
        if abs(d - 1.0) > tol:
            bad.append((i, j, d))
    return bad

def write_xyz(data, filename):
    N = data.shape[0]
    with open(filename, 'w') as f:
        f.write(str(N) + '\n\n')
        for i in range(N):
            f.write(f'1\t{data[i, 0]}\t{data[i, 1]}\t{data[i, 2]}\n')

def generate_chain_with_fixed_distance(point1, point2, num_points, tol=1e-3, max_iter=50):
    """Uses binary search to find a Bezier control point that yields a specific arc length."""
    p1 = np.array(point1, dtype=float)
    p2 = np.array(point2, dtype=float)

    if num_points < 2:
        raise ValueError("num_points must be >= 2")

    target_length = (num_points - 1) * 1.0
    vec = p2 - p1
    dist = np.linalg.norm(vec)
    if dist == 0:
        return np.tile(p1, (num_points, 1))

    dir_unit = vec / dist
    midpoint = (p1 + p2) / 2

    low, high = 0.0, 5 * target_length 
    best_control_offset = 0.0

    for _ in range(max_iter):
        mid = (low + high) / 2
        z_axis = np.array([0, 0, 1])
        perp = np.cross(dir_unit, z_axis)
        if np.linalg.norm(perp) < 1e-6:
            x_axis = np.array([1, 0, 0])
            perp = np.cross(dir_unit, x_axis)
        perp_unit = perp / np.linalg.norm(perp)

        control_point = midpoint + perp_unit * mid
        length, _, _ = bezier_length(p1, control_point, p2)
        diff = length - target_length

        if abs(diff) < tol:
            best_control_offset = mid
            break

        if diff > 0:
            high = mid
        else:
            low = mid
        best_control_offset = mid

    control_point = midpoint + perp_unit * best_control_offset
    points = sample_points_by_arc_length(p1, control_point, p2, num_points)
    return points[1:-1]

def rotate_to_principal_axis(knotcore: np.ndarray):
    """Rotates the coordinates based on the principal moments of inertia."""
    X = knotcore
    X -= rmsd.centroid(X)

    I = np.zeros((3, 3))
    for x in X:
        I[0, 0] += x[1]**2 + x[2]**2
        I[1, 1] += x[0]**2 + x[2]**2
        I[2, 2] += x[0]**2 + x[1]**2
        I[0, 1] -= x[0] * x[1]
        I[0, 2] -= x[0] * x[2]
        I[1, 2] -= x[1] * x[2]

    I[1, 0] = I[0, 1]
    I[2, 0] = I[0, 2]
    I[2, 1] = I[1, 2]
    
    eigvals, eigvecs = np.linalg.eig(I)
    idx = eigvals.argsort()
    eigvecs = eigvecs[:, idx]

    if np.linalg.det(eigvecs) < 0:
        eigvecs[:, 2] = -eigvecs[:, 2]

    X_copy = np.dot(X, eigvecs)
    return X_copy

def write_lammps(data, filename, chain_type="open", Lx=200, Ly=200, Lz=200):
    """Writes system data in LAMMPS format."""
    N = data.shape[0]
    with open(filename, 'w') as f:
        f.write("# LAMMPS input file\n")
        f.write(f'{N} atoms\n')
        
        # Bond count
        f.write(f'{N if chain_type=="close" else N-1} bonds\n')
        # Angle count
        f.write(f'{N if chain_type=="close" else N-2} angles\n')

        f.write('\n1 atom types\n1 bond types\n1 angle types\n')

        # Box dimensions
        f.write(f'\n{-Lx} {max(Lx, np.max(data[:,0]))} xlo xhi\n')
        f.write(f'{-Ly} {max(Ly, np.max(data[:,1]))} ylo yhi\n')
        f.write(f'{-Lz} {max(Lz, np.max(data[:,2]))} zlo zhi\n')

        f.write('\nMasses\n\n1 1.0\n')

        f.write('\nAtoms\n\n')
        for i in range(N):
            f.write(f'{i+1}\t1\t1\t{data[i,0]}\t{data[i,1]}\t{data[i,2]}\n')

        f.write('\nBonds\n\n')
        for i in range(N-1):
            f.write(f'{i+1}\t1\t{i+1}\t{i+2}\n')
        if chain_type == "close":
            f.write(f'{N}\t1\t{N}\t1\n')

        f.write('\nAngles\n\n')
        for i in range(N-2):
            f.write(f'{i+1}\t1\t{i+1}\t{i+2}\t{i+3}\n')
        if chain_type == "close":
            f.write(f'{N-1}\t1\t{N-1}\t{N}\t1\n')
            f.write(f'{N}\t1\t{N}\t1\t2\n')

def bezier_generator_700(core, guide, total):
    """Connects a core structure and guide structure using two long Bezier bridges."""
    core = core - 60
    len1, len2 = len(core), len(guide)
    n = total - len1 - len2
    points1_num = 350
    points2_num = n - points1_num
    
    points1 = generate_chain_with_fixed_distance(guide[-1], core[0], num_points=points1_num + 2)
    points2 = generate_chain_with_fixed_distance(core[-1], guide[0], num_points=points2_num + 2)
    
    return np.vstack([guide, points1, core, points2])

def bezier_generator(core, guide, total):
    """Connects a core structure and guide structure using two Bezier bridges."""
    core = core - 50
    len1, len2 = len(core), len(guide)
    n = total - len1 - len2
    points1_num = 100
    points2_num = n - points1_num
    
    points1 = generate_chain_with_fixed_distance(guide[-1], core[0], num_points=points1_num + 2)
    points2 = generate_chain_with_fixed_distance(core[-1], guide[0], num_points=points2_num + 2)
    
    return np.vstack([guide, points1, core, points2])

# Execution Loop
for knot in [0, 31]:
    for length in [300]:
        core_coords = np.loadtxt(f"core/pos_knot{knot}.txt")
        guide_coords = np.loadtxt("curve_150.txt")
        
        final_data = bezier_generator(core_coords, guide_coords, length)
        
        write_xyz(final_data, f'MC_{knot}_L{length}_close_rigid.xyz')
        write_lammps(final_data, f"lammps_{knot}_L{length}_close.data", 
                    chain_type="close", Lx=length, Ly=length, Lz=length)
        
        errors = check_distance(final_data)
        if errors:
            print(f"Found {len(errors)} distance anomalies:")
            for i, j, d in errors:
                print(f"  Atoms {i} and {j}, distance = {d:.4f}")
        else:
            print(f"Knot {knot}, Length {length}: All distances valid.")
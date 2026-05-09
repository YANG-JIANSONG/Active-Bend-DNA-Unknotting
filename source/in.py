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
        f"read_data lammps_1006_{type}_{L}_{degree}.data",
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
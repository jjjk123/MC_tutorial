from typing import List, Set, Dict, Tuple, Optional
import random
import math
import numpy as np

def generate(N, d_min:float, L:float) -> np.ndarray:
    """ Initializes positions of all atoms in 2D; points (x, y) should not overlap
    both x and y in the range [0,L]
    """
    n = int(math.sqrt(N))
    x_list = np.arange(n, L, L/n)
    y_list = np.arange(n, L, L/n)
    atoms = np.array(np.meshgrid(x_list, y_list)).T.reshape(-1, 2)
    return atoms

def print_pdb(atoms, model) ->None:
    """ Prints all positions on a screen"""
    pdb = ""
    pdb += "MODEL    " + str(model) + '\n'
    for index in range(len(atoms)):
        pdb += str("ATOM   %4d%4s  AAA A%4d    %8.3f%8.3f%8.3f  0.50 35.88           A" % \
            (index, "AR", index, atoms[index][0], atoms[index][1], 0)) + '\n'
    return pdb


def energy(atoms) -> float:
    """
    Evaluates energy, e.g. square well
    """
    en = 0
    for i in range(len(atoms)):
        en += energy(atoms, i, 1, 4)
    return en


def energy(atoms, index:int, d_min:float, d_max:float) -> float:
    """
    Evaluates energy of a single atom
    """
    for i in range(len(atoms)):
        d = math.sqrt((atoms[index-1][0] - atoms[index][0])**2 + (atoms[index][1] - atoms[index][1])**2)
    en = 0
    if d < d_min:
        return 10000000000
    elif d < d_max:
        en += -1
    else:
        return 0
    return en
    # PWB


def metropolis(atoms, T:float) -> bool:
    """
    1) Randomly selects an atom (i_moved)
    6) e_before = energy(atoms, i_moved)
    2) Randoms dx, dy, say [-0.5A, 0.5A]
    3) copy position of i_moved
    4) atoms[i_moved][0] += dx
       atoms[i_moved][1] += dy
    5) Apply PBC
    7) e_after = energy(atoms, i_moved)
    8) Metropolis criterion
    """
    i_moved = int(random.random()*10)
    en_before = energy(atoms, i_moved)
    dx = random.uniform(-0.5, 0.5)
    dy = random.uniform(-0.5, 0.5)
    prev_x, prev_y = atoms[i_moved]
    atoms[i_moved][0] += dx
    atoms[i_moved][1] += dy

    if atoms[i_moved][0] > L:
        atoms[i_moved][0] = L - dx
    elif atoms[i_moved][1] > L:
        atoms[i_moved][1] = L - dy

    en_after = energy(atoms, i_moved)

    if en_after < en_before:
        pass
    else:
        r = random.random()
        if math.exp(-(en_after - en_before)/T) < r:
            pass
        else:
            atoms[i_moved][0] = prev_x
            atoms[i_moved][1] = prev_y

if __name__ == "__main__":

    N_atoms: int = 9           # --- the number of atoms (particles)
    L = 10                      # --- size of the periodic box in A
    r0: float = 2.0             # --- atom radius
    T:float = 1.0               # --- temperature of the simulation

    # --- positions
    atoms: List[List] = [[0, 0] for i in range(N_atoms)]
    atoms = generate(N_atoms, r0, L)
    print(atoms)
    pdb = print_pdb(atoms, 1)
    print(pdb)
    f = open('MC_simulation.pdb', 'w')
    f.write(pdb)
    for i in range(100):
        metropolis(atoms,T)
        print_pdb(atoms)

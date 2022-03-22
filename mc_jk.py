from typing import List, Set, Dict, Tuple, Optional
import random
import math
import numpy as np

def generate(N, L:float) -> np.ndarray:
    """ Initializes positions of all atoms in 2D; points (x, y) should not overlap
    both x and y in the range [0,L]
    """
    n = int(math.sqrt(N))
    x_list = np.arange(n, L, L/n)
    y_list = np.arange(n, L, L/n)
    atoms = np.array(np.meshgrid(x_list, y_list)).T.reshape(-1, 2)
    return atoms

def print_pdb(atoms, model, file_name) ->None:
    """ Prints all positions on a screen"""
    pdb = ""
    pdb += "MODEL    " + str(model) + '\n'
    for index in range(len(atoms)):
        pdb += str("ATOM   %4d%4s  AAA A%4d    %8.3f%8.3f%8.3f  0.50 35.88           A" % \
            (index, "AR", index, atoms[index][0], atoms[index][1], 0)) + '\n'
    pdb += "ENDMDL\n"
    with open(file_name, 'a') as f:
        f.write(pdb)

def total_energy(atoms, L, d_min, d_max) -> float:
    """
    Evaluates energy, e.g. square well
    """
    en = 0
    for i in range(len(atoms)):
        en += energy(atoms, L, i, d_min, d_max)
    return en


def energy(atoms, L, index:int, d_min:float, d_max:float) -> float:
    """
    Evaluates energy of a single atom
    """
    en = 0
    for i in range(len(atoms)):
        if i == index:
            continue
        dx = abs(atoms[index][0] - atoms[i][0])
        dy = abs(atoms[index][1] - atoms[i][1])
        dx = dx if dx < L/2 else L - dx 
        dy = dy if dy < L/2 else L - dy
        d = math.sqrt(dx**2 + dy**2)
        # if index == 0:
        #     print(i, dx, dy, d)
        if d < d_min:
            en += 10000000000
        elif d < d_max:
            en += -1
    return en

def metropolis(atoms, L, T:float, d_min, d_max) -> bool:
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
    N = len(atoms)
    i_moved = int(random.random() * N)
    en_before = energy(atoms, L, i_moved, d_min, d_max)
    dx = random.uniform(-0.5, 0.5)
    dy = random.uniform(-0.5, 0.5)
    prev_x, prev_y = atoms[i_moved]
    atoms[i_moved][0] += dx
    atoms[i_moved][1] += dy

    if atoms[i_moved][0] > L:
        atoms[i_moved][0] = L - atoms[i_moved][0]
    elif atoms[i_moved][1] > L:
        atoms[i_moved][1] = L - atoms[i_moved][1]

    en_after = energy(atoms, L, i_moved, d_min, d_max)

    d_en = en_after - en_before 
    if d_en < 0:
        pass
    else:
        r = random.random()
        if math.exp(-(d_en)/T) < r:
            pass
        else:
            atoms[i_moved][0] = prev_x
            atoms[i_moved][1] = prev_y

if __name__ == "__main__":

    N_atoms: int = 100           # --- the number of atoms (particles)
    L = 30                     # --- size of the periodic box in A
    T:float = 1.0               # --- temperature of the simulation
    d_min = 2
    d_max = 3

    # --- positions
    atoms = generate(N_atoms, L)
    en_prev = total_energy(atoms, L, d_min, d_max)
    print('prev_en', en_prev)
    for j in range(100):
        for i in range(100):
            metropolis(atoms, L, T, d_min, d_max)
        print(total_energy(atoms, L, d_min, d_max))
        print_pdb(atoms, i, 'MC_simulation.pdb')
    en_after = total_energy(atoms, L, d_min, d_max)
    print('after_en', en_after)



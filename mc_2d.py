from typing import List, Set, Dict, Tuple, Optional

def generate(atoms, d_min:float, L:float) -> None:
    """ Initializes positions of all atoms in 2D; points (x, y) should not overlap
    both x and y in the range [0,L]
    """
    pass


def print_pdb(atoms) ->None:
    """ Prints all positions on a screen"""
    pass


def energy(atoms) -> float:
    """
    Evaluates energy, e.g. square well
    """
    pass


def energy(atoms, index:int) -> float:
    """
    Evaluates energy of a single atom
    """
    pass


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


if __name__ == "__main__":

    N_atoms: int = 10           # --- the number of atoms (particles)
    L = 20                      # --- size of the periodic box in A
    r0: float = 2.0             # --- atom radius
    T:float = 1.0               # --- temperature of the simulation

    # --- positions
    atoms: List[Tuple[int, int]] = [(0, 0) for i in range(N_atoms)]
    generate(atoms, r0, L)
    for i in range(100):
        metropolis(atoms,T)
        print_pdb(atoms)

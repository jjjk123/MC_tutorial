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

atoms = [[3.,         3.        ],
[3.         ,6.33333333],
 [3.         ,9.66666667],
 [6.33333333 ,3.        ],
 [6.33333333 ,6.33333333],
 [6.33333333 ,9.66666667],
 [9.66666667 ,3.        ],
 [9.66666667 ,6.33333333],
 [9.66666667 ,9.66666667]]

 print(atoms)
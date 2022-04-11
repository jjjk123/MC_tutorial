# x = -0.5
# y = 0.5

# if x < 0:
#     print('hello')
# if y < 0:
#     print('hello')
# else:
#     print('a')


# l = [1, 2, 3]
# for i, x in enumerate(l):
#     print(i, x)

import numpy as np

with open('MC_simulation.pdb', 'r') as f:
    file = f.read()

atoms = []
last_model = file[file.rfind('MODEL'):].splitlines()
for atom in last_model[1:-1]:
    x = float(atom[32:39].replace(' ', ''))
    y = float(atom[40:47].replace(' ', ''))
    atoms.append([x, y])
print(atoms)
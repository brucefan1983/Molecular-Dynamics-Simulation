from ase import Atoms
from ase.calculators.tip4p import angleHOH, rOH
from ase.io import write
import numpy as np
import scipy.constants as spc

# Set up water box at 20 deg C density
x = angleHOH * np.pi / 180 / 2
positions = [[0, 0, 0],
                 [0, rOH * np.cos(x), rOH * np.sin(x)],
                 [0, rOH * np.cos(x), -rOH * np.sin(x)]]
atoms = Atoms('OH2', positions=positions)

# Denote number of molecules, molar mass and target density
number_of_molecules = 512
molar_mass_water = 18.01528
target_density = 0.9994
NA = spc.value('Avogadro constant')
box_length = ((molar_mass_water / NA) / (target_density / 1e24))**(1 / 3.)
atoms.set_cell((box_length, box_length, box_length))
atoms.center()

number_or_replica = int(np.ceil(number_of_molecules ** (1 / 3.)))
atoms = atoms.repeat(number_or_replica)
atoms.set_pbc(True)

# Write water box
write('model.xyz', atoms)

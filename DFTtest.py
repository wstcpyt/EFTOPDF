__author__ = 'yutongpang'
from PyQuante.Molecule import Molecule
h2 = Molecule('H2',
                 [(1,  (0.00000000,     0.00000000,     0.36628549)),
                  (1,  (0.00000000,     0.00000000,    -0.36628549))],
                 units='Angstrom')
from PyQuante.dft import dft
en,orbe,orbs = dft(h2)
print "HF Energy = ",orbs
"""
Created on 2023-01-17 16:12:01.409428
@author: suvobanik
"""

import math
import os

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.io.cif import CifWriter
from pymatgen.io.vasp import Poscar
from tqdm import tqdm

# In[10]:


class StructureWriter(object):
    def __init__(
        self,
        dumpfile,
        outpath="structures",
        objfile="energy.dat",
        file_format="poscar",
    ):

        """
        :dumpfile: full path of the file, where structures are written by MCTS.
        :outpath: full path of the directory where the structures are to be extracted.
        :objfile: full path of the file, where the objectives are to be written.
        :file_format: file format of the extracted structures. Currently only cif or poscar is supported .
        """

        self.dumpfile = dumpfile
        self.outpath = outpath
        self.objfile = objfile
        if file_format == "poscar":
            self.writer = Poscar
        elif file_format == "cif":
            self.writer = CifWriter
        else:
            raise Exception("The file format not specified")

        if not os.path.exists(self.outpath):
            os.mkdir(self.outpath)

        try:
            os.remove(self.objfile)
        except:
            pass

    def gel_latt_coords(self, parameters):
        lattice = parameters[:6]
        coords = parameters[6:]
        coords = coords.reshape(int(coords.shape[0] / 3), 3)
        return lattice, coords

    def write(self, num_to_write, sort=True):
        structurelist = []

        with open(self.dumpfile, "r") as infile:
            for i, line in enumerate(infile):
                col = line.split("|")
                energy = float(col[2])
                if math.isnan(energy):
                    continue
                parms = np.array([float(val) for val in col[0].split()])
                species = col[1].split()
                structurelist.append([energy, parms, species])

        if sort:
            structurelist.sort(key=lambda x: x[0])

        for i, val in tqdm(enumerate(structurelist[:num_to_write])):

            parameters = val[1]
            species = val[2]
            lattice, coords = self.gel_latt_coords(parameters)

            lattice = Lattice.from_parameters(
                a=lattice[0],
                b=lattice[1],
                c=lattice[2],
                alpha=lattice[3],
                beta=lattice[4],
                gamma=lattice[5],
            )

            struct = Structure(lattice, species, coords, to_unit_cell=True)

            with open(self.objfile, "a") as outfile:
                outfile.write("{} {}\n".format(i, val[0]))

            Poscar(struct).write_file("{}/{}.POSCAR".format(self.outpath, i))

import os
import sys
import numpy as np
from pymatgen.io.vasp import Poscar
from pymatgen.core import Structure,Lattice
import math


path = sys.argv[1]
number = int(sys.argv[2])


if not os.path.exists(path):
    os.mkdir(path)




structurelist = []


def gel_latt_coords(parameters):
    lattice = parameters[:6]
    coords = parameters[6:]
    coords = coords.reshape(int(coords.shape[0] / 3), 3)
    return lattice, coords



with open("dumpfile.dat","r") as infile:
    for i,line in enumerate(infile):
        col = line.split("|")

        parms = np.array([float(val) for val in col[0].split()])
        species = col[1].split()
        energy = float(col[2])
        
        if math.isnan(energy):
            continue


        structurelist.append([energy,parms,species])


structurelist.sort(key=lambda x:x[0])




for i,val in enumerate(structurelist[:number]):
    
    parameters = val[1]
    species = val[2]
    lattice,coords = gel_latt_coords(parameters)
    
    
    lattice = Lattice.from_parameters(a=lattice[0],b=lattice[1],
                                      c=lattice[2],alpha=lattice[3],
                                      beta=lattice[4],gamma=lattice[5])
    
    
    
    struct = Structure(lattice, species ,coords,to_unit_cell=True)
    

    with open("energy.dat","a") as outfile:
        outfile.write("{} {}\n".format(i,val[0]))

    Poscar(struct).write_file("{}/{}.POSCAR".format(path,i))

    print("structure {} created in {}".format(i,path))





















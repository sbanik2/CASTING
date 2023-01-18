from CASTING.writer import StructureWriter


num_to_write = 10 # number of stuctures to extract
writer = StructureWriter("dumpfile.dat",outpath="structures",objfile="energy.dat",file_format="poscar") # "poscar" or "cif"
writer.write( num_to_write, sort=True) # sort to arrange in increasing order of energy







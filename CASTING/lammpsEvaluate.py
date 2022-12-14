#!/usr/bin/env python
# coding: utf-8

# In[1]:

"""
Created on 2022-12-13 20:36:01.409428
@author: suvobanik
"""


import math

from lammps import lammps
from pymatgen.io.lammps.data import LammpsData

from CASTING.clusterfun import check_constrains, parm2struc, struc2param

# In[ ]:

cmds = ["-screen", "log.screen"]
lmp = lammps(cmdargs=cmds)
# lmp  = lammps()


class LammpsEvaluator(object):
    """
    Performs energy evaluations on a offspring
    """

    def __init__(
        self,
        constrains,
        pair_style,
        pair_coeff,
    ):

        self.constrains = constrains
        self.pair_style = pair_style
        self.pair_coeff = pair_coeff

    def evaluate(self, structData):

        if not check_constrains(structData, self.constrains, verbose=False):
            return structData, 1e300

        struct = parm2struc(structData, self.constrains)
        LammpsData.from_structure(struct, atom_style="atomic").write_file("in.data")

        lmp.command("clear")
        lmp.command("dimension 3")
        lmp.command("box tilt large")
        lmp.command("units metal")
        lmp.command("atom_style atomic")
        lmp.command("neighbor 2.0 bin")
        lmp.command("atom_modify map array sort 0 0")
        lmp.command("boundary f f f")
        lmp.command("read_data in.data")
        lmp.command("{}".format(self.pair_style))
        lmp.command("{}".format(self.pair_coeff))
        lmp.command("thermo 1000")
        lmp.command("thermo_style custom step etotal atoms vol")
        lmp.command("thermo_modify format float %5.14g")
        lmp.command("variable potential equal pe/atoms")
        lmp.command("neigh_modify one 5000 delay 0 every 1 check yes")
        lmp.command("run 0 pre no")

        # ------------guard for bad structures---------

        energy = lmp.extract_variable("potential", None, 0)

        if math.isinf(float(energy)):
            return structData, 1e300
        elif math.isnan(float(energy)):
            return structData, 1e300
        else:
            pass
        # ---------------------------------------------

        lmp.command("minimize 1.0e-8 1.0e-8 10000 10000")
        lmp.command("write_data min.geo")
        lmp.command("run 0 pre no")

        energy = lmp.extract_variable("potential", None, 0)

        minstruct = LammpsData.from_file("min.geo", atom_style="atomic").structure

        minData, mineng = struc2param(
            minstruct,
            energy,
            self.constrains,
            CheckFrConstrains=True,
            writefile="dumpfile.dat",
        )

        return minData, mineng

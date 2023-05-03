from ase import Atoms
from ase.calculators.eam import EAM
from ase.optimize import BFGS
from mendeleev import element
from pymatgen.core import Lattice, Structure
from pymatgen.core.periodic_table import Element
from pymatgen.io.ase import AseAtomsAdaptor

from CASTING.clusterfun import check_constrains, parm2struc, struc2param


def structure_to_cell(structure):
    lattice = structure.lattice.matrix.tolist()
    positions = [list(site.coords) for site in structure.sites]
    numbers = [Element(site.specie.symbol).Z for site in structure.sites]
    return (lattice, positions, numbers)


def cell_to_structure(cell):
    (lattice, positions, numbers) = cell
    return Structure(
        Lattice(lattice), [element(int(num)).symbol for num in numbers], positions
    )


def structure_to_Ase(structure):
    (lattice, positions, numbers) = structure_to_cell(structure)

    return Atoms(numbers, positions=positions, cell=lattice, pbc=[0, 0, 0])


class AseEamEvaluator:
    def __init__(
        self,
        constrains,
        par_file,
    ):

        self.constrains = constrains
        self.par_file = par_file

    def evaluate(self, structData):

        if not check_constrains(structData, self.constrains, verbose=False):
            return structData, 1e300

        struct = parm2struc(structData, self.constrains)

        calc = EAM(
            potential=self.par_file,
            elements=["Au"],
            Z=[79],
            lattice=["fcc"],
            mass=[196.966570],
            a=[11.999000],
            form="eam",
        )

        atoms = structure_to_Ase(struct)
        atoms.calc = calc
        mininimizer = BFGS(atoms, maxstep=0.2, logfile=None)
        mininimizer.run(fmax=0.001, steps=200)
        atoms.set_constraint()

        minstruct = AseAtomsAdaptor.get_structure(atoms)

        minData, mineng = struc2param(
            minstruct,
            atoms.get_potential_energy(),
            self.constrains,
            CheckFrConstrains=True,
            writefile="dumpfile.dat",
        )

        return minData, mineng

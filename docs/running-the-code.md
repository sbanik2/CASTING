### Running the code
<p align="justify"> First, all the parameters crystal (constrains), LAMMPS parameters (pair style, pair coefficient etc.) and the perturbation parameter need to be set.  The composition is given for e.g., a Au<2>Al<3> as "composition":{"Au":2,"Al:3"}. In a file (for e.g., RunOpt.py) we define, </p>


``` python
from CASTING.utilis import r_datafame,get_lattice
from CASTING.MCTS import MCTS
from CASTING.clusterfun import createRandomData
from CASTING.lammpsEvaluate import LammpsEvaluator
from CASTING.perturb import perturbate


# In[6]:

import random
import numpy as np

seed = 12

random.seed(seed)
np.random.seed(seed)



# In[9]:


r_min = {"Au-Au":2}       # minimum allowed interatomic distance
r_max = {"Au-Au":4}      # maximum allowed interatomic distance
box_dim = 50                  # box dimension 

#--------crystal constrains-----------

constrains = {
    "composition":{"Au":1},
    "atoms":13,
    "r_min":r_datafame(r_min),
    "r_max":r_datafame(r_max),
    "lattice":get_lattice(box_dim),  
    }




#-------------perturbation------------


pt = {
    'max_mutation': 0.05,  # mutation  in fraction of  box length
}



#-----------------lammps parameters-------------------

lammps_par = {
    'constrains':constrains,
    'pair_style': "pair_style eam",
    'pair_coeff': "pair_coeff * * Au.eam"   # Provide full path of potential file here
}
    
```

Once the parameters are set, the evaluator (LAMMPS calculator) & the perturbator is to be initialized and a structure for root node is created.

 ``` python

rootdata = createRandomData(constrains,multiplier= 10)

perturb = perturbate(**pt). perturb
evaluator = LammpsEvaluator(**lammps_par).evaluate 

```

Finally, Call MCTS with all the hyperparameters added. Details of individual hyperparameters for the optimizer can be found here [Paper](https://doi.org/10.48550/arXiv.2212.12106).
        
 ``` python
        
MCTS(
    rootdata,
    perturb ,
    evaluator,
    niterations=2000,
    headexpand=10,
    nexpand=3,
    nsimulate=3,
    nplayouts=10,
    exploreconstant=1,
    maxdepth=12,
    a=0,
    selected_node=0,
)
        
```

Run the code 

```
python RunOpt.py 
```

The optimization produces a "dumpfile.dat" output containing all the crystal parameters and the energy values as the output. To extract the structures in either 'poscar' or 'cif' format, one can use the  'StructureWriter' module.

 ``` python
 
from CASTING.writer import StructureWriter

num_to_write = 10 # number of stuctures to extract
writer = StructureWriter(
                         "dumpfile.dat",
                          outpath="structures",
                          objfile="energy.dat",
                          file_format="poscar" # "poscar" or "cif"
                          ) 
writer.write( num_to_write, sort=True) # sort to arrange in increasing order of energy
 
 
 ```
This will extract 'num_to_write' number of structures in ascending order of objective. 


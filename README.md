#CASTING 

<a name="readme-top"></a>



[![Release][release-shield]]
[![License][license-shield]][license-url]
[![Downloads][download-shield]]
![Commit][commit-shield]
[![Size][size-shield]]
[![DOI][DOI-shield]][DOI-url]



# 

<p align="justify"> A Continuous Action Space Tree search for INverse desiGn (CASTING) Framework and Materials Discovery</p>


## Table of Contents
- [Introduction](#Introduction)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Running the code](#running-the-code)
- [Optimization of Gold nanocluster](#optimization-of-Gold-nanocluster)
- [Carbon metastable polymorphs](#carbon-metastable-polymorphs)
- [Citation](#citation)
- [License](#license)


<p align="right">(<a href="#readme-top">back to top</a>)</p>


## Introduction

<p align="justify">A pseudocode implementation of CASTING framework (https://doi.org/10.48550/arXiv.2212.12106) for optimization of atomic nanoclusters only. This code uses MCTS (Monte Carlo Tree Search) as base optimizer.<br/> 
&emsp;&emsp;&emsp;&emsp;Fast and accurate prediction of optimal crystal structure, topology, and microstructures is important for accelerating the design and discovery of new materials. Material properties are strongly correlated to the underlying structure and topology – inverse design is emerging as a powerful tool to discover new and increasingly complex materials that meet targeted functionalities. CASTING provides a unified framework for fast, scalable and accurate prediction & design of materials. </p>


<p align="center"> <a href="url"><img src="https://github.com/sbanik2/CASTING/blob/main/figs/MCTS.png?raw=true" align="center" height="500" width="600" ></a> </p>


<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Prerequisites
This package requires:
- [scipy](https://scipy.org/)
- [LAMMPS](https://www.lammps.org/)
- [pymatgen](https://pymatgen.org/)
- [pandas](https://pandas.pydata.org/)
- [numpy](https://numpy.org/)
- [networkx](https://networkx.org/)
- [ase](https://wiki.fysik.dtu.dk/ase/#)

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Installation

### Manual Installation
[Install the anaconda package](https://docs.anaconda.com/anaconda/install/). Then, 
```
conda env create --name CASTING
conda activate CASTING
git clone https://github.com/sbanik2/CASTING.git
pip install -r requirement.txt
python setup.py install
```

### Installation with pypi

```
pip install CASTING

```
***The package requires python lammps binding to run. First, lammps package needs to be downloaded from [LAMMPS download] (https://www.lammps.org/download.html) and compiled. The instructions on python integration can be found here [LAMMPS-Python] (https://docs.lammps.org/Python_install.html).

<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Running the code
<p align="justify"> First all the parameters crystal (constrains), LAMMPS parameters (pair style, pair coefficient etc.) and the perturbation parameter need to be set.  The composition is given for e.g., a Au<2>Al<3> as "composition":{"Au":2,"Al:3"}. In a file (for e.g., RunOpt.py) we define, </p>


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

Finally, Call MCTS with all the hyperparameters added. Details of individual hyperparameters for the optimizer can be found here (
https://doi.org/10.48550/arXiv.2212.12106).
        
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

The optimization produces a "dumpfile.dat" output containing all the crystal parameters and the energy values as the output.  A script “Createstruct.py” for extracting the structure in POSCAR format is given in the “post” directory. To run this

```
python Createstruct.py <path-to-extraction-directory> <number-of-structure-to-extract>
```
This will extract <number-of-structure-to-extract> number of structures in ascending order of objective.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Optimization of Gold nanocluster

An example optimization of Gold (Au) nanocluster is given in "example" directory. We have used CASTING to optimize the already known global minima (Sutton-Chen) of 13 atom Au nanocluster (Icosahedral structure).  Details on additional examples can be found in (
https://doi.org/10.48550/arXiv.2212.12106)

<p align="center"> <a href="url"><img src="https://github.com/sbanik2/CASTING/blob/main/figs/sutton_chen.gif" align="center" height="200" width="200" ></a> </p>

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Carbon metastable polymorphs

We have also used CASTING to sample metastable polymorphs of Carbon(C). All the structures are then further relaxed with DFT. The unique polymorphs and their corresponding DFT energies have been provided in “C_polymorphs” directory.

<p align="center"> <a href="url"><img src="https://github.com/sbanik2/CASTING/blob/main/figs/MetastableC.png" align="center" height="400" width="500" ></a> </p>

## Citation
```
Banik, S., Loefller, T., Manna, S., Srinivasan, S., Darancet, P., Chan, H., ... & Sankaranarayanan, S. K. (2022). A Continuous Action Space Tree search for INverse desiGn (CASTING) Framework for Materials Discovery. arXiv preprint arXiv:2212.12106.

```
    
<p align="right">(<a href="#readme-top">back to top</a>)</p>
        
## License
CASTING is distributed under MIT License. See `LICENSE` for details.
    
    
<p align="right">(<a href="#readme-top">back to top</a>)</p>  
    
<!--LINKS -->


[release-shield]:https://img.shields.io/github/v/release/sbanik2/CASTING
[license-shield]: https://img.shields.io/github/license/sbanik2/CASTING
[license-url]: https://github.com/sbanik2/CASTING/blob/main/LICENSE
[download-shield]: https://img.shields.io/github/downloads/sbanik2/CASTING/total
[commit-shield]::https://img.shields.io/github/last-commit/sbanik2/CASTING
[size-shield] : https://img.shields.io/github/languages/code-size/sbanik2/CASTING
[DOI-shield]: https://img.shields.io/badge/DOI-doi.org/10.48550/arXiv.2212.12106-red
[DOI-url]:https://doi.org/10.48550/arXiv.2212.12106
    
    

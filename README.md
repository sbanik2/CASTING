# ClusterOpt

<p align="justify"> Implementation of Basin Hopping Global Optimization method for optimization of Atomic nanoclusters. </p>


## Table of Contents
- [Introduction](#Introduction)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Running the code](#Running-the-code)
- [Citation](#data-availability)
- [License](#license)

## Introduction
This package implements the Basin Hopping global optimization method in [scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.basinhopping.html) along with [LAMMPS](https://www.lammps.org/) MD (Molecular dynamic simulation) package to find the global minima of atomic nanoclusters (Single or Multicomponent)

<p align="center"> <a href="url"><img src="https://github.com/sbanik2/ClusterOpt/blob/main/Figs/cluster_opt.png" align="center" height="300" width="300" ></a> </p>



## Prerequisites
This package requires:
- [scipy](https://scipy.org/)
- [LAMMPS](https://www.lammps.org/)
- [pymatgen](https://pymatgen.org/)
- [pandas](https://pandas.pydata.org/)
- [numpy](https://numpy.org/)


## Installation

### Manual Installation
[Install the anaconda package](https://docs.anaconda.com/anaconda/install/). Then, 
```
conda env create --name ClusterOpt -f environment.yml
conda activate ClusterOpt
git clone https://github.com/sbanik2/ClusterOpt.git
python setup.py install
```

### Installation with pypi

```
pip install ClusterOpt

```
***The package requires python lammps binding to run. First, lammps package needs to be downloaded from  [LAMMPS download](https://www.lammps.org/download.html) and compiled. The instructions on python integration can be found here [LAMMPS-Python](https://docs.lammps.org/Python_install.html).


### Running the code
<p align="justify"> An example of running run directory provided in the example section. First all the parameters crystal and the lammps pair_style and pair_coeff should be set. The composition is given for e.g., a Au<2>Al<3> as "composition":{"Au":2,"Al:3"}, the minimum interatomic distances as a pandas data frame with rows and columns belonging to each species in the same order they are mentioned in the composition. E.g.,</p>

``` python
import numpy as np
import pandas as pd
from ClusterOpt.Evaluator import Custom_minimize, LammpsEvaluator
from ClusterOpt.Structfunc import createRandomData
from ClusterOpt.utilis import CreateStructure,Status
from scipy.optimize import basinhopping
from scipy.optimize import minimize


# minimum distance criteria between the atoms

r = pd.DataFrame(np.array([[1.3]]),columns=["Au"],index= ["Au"])

constrains = {     
        "composition":{"Au":1},
        "atoms":13,
        "vpa":[16, 20],
        "r":r,       
        }


lammps_args = {
        "pair_style" : "pair_style eam",
        "pair_coeff" : "pair_coeff * * Au.eam",
        "pad":20
            }

args = (constrains,lammps_args)

```

Once the parameters are set, the LammpsEvaluator can be initialized and the basinhopping optimizer can be set for the run

 ``` python

structData = createRandomData(constrains,trials = 1000)
x0 = structData["parameters"]


fun = LammpsEvaluator(constrains,lammps_args).snapshot
res = minimize(fun, x0, args=args, method=Custom_minimize)

print(res.x,res.fun)

minimizer_kwargs = {"method":Custom_minimize,"args":args}

def print_fun(x, f, accepted):
    print("at minimum %.2e accepted %d" % (f, int(accepted)))

basinhopping(fun,
             x0, 
             niter=1000,
             T=1.0, 
             stepsize=0.5, 
             minimizer_kwargs=minimizer_kwargs, 
             take_step=None, 
             accept_test=None, 
             callback=print_fun, 
             interval=50, 
             disp=False, 
             niter_success=None, 
             seed=None, 
            )
```
details of individual hyperparameters for the basin hopping optimizer can be found here [BasinHopping](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.basinhopping.html).
The optimization produces a "dumpfile.dat" output containing all the crystal parameters and the energy values as the output. The best values and the best crystal structure corresponding to the least energy can be extracted using
        
 ``` python
        
Status("dumpfile.dat")

CreateStructure("dumpfile.dat","./",nStructure=1)
        
```


        
### Referenece
   
1. Wales, D. J., & Doye, J. P. (1997). Global optimization by basin-hopping and the lowest energy structures of Lennard-Jones clusters containing up to 110 atoms. The Journal of Physical Chemistry A, 101(28), 5111-5116.
2. Wales, D. J., & Scheraga, H. A. (1999). Global optimization of clusters, crystals, and biomolecules. Science, 285(5432), 1368-1372.   


### Citation
```
@article{manna2022learning,
  title={Learning in continuous action space for developing high dimensional potential energy models},
  author={Manna, Sukriti and Loeffler, Troy D and Batra, Rohit and Banik, Suvo and Chan, Henry and Varughese, Bilvin and Sasikumar, Kiran and Sternberg, Michael and Peterka, Tom and Cherukara, Mathew J and others},
  journal={Nature communications},
  volume={13},
  number={1},
  pages={1--10},
  year={2022},
  publisher={Nature Publishing Group}
}
```
        
### License
ClusterOpt is licensed under the MIT License


#!/usr/bin/env python
# coding: utf-8

# In[1]:



# In[2]:


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


r_min = {"Au-Au":2}
r_max = {"Au-Au":4}
box_dim = 50 # Make sure there is enough vaccum in the sides

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
    'max_mutation': 0.05,  #Put in fraction of the box length 0.01 means 100*0.01 =1Angs
}



#-----------------lammps parameters-------------------

lammps_par = {
    'constrains':constrains,
    'pair_style': "pair_style eam",
    'pair_coeff': "pair_coeff * * Au.eam"   
}
    


    


# In[10]:


rootdata = createRandomData(constrains,multiplier= 10)

perturb = perturbate(**pt). perturb
evaluator = LammpsEvaluator(**lammps_par).evaluate


# In[ ]:


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



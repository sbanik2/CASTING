"""
Created on 2022-12-13 20:36:01.409428
@author: suvobanik
"""


from random import random

import numpy as np


class perturbate(object):
    def __init__(
        self,
        max_mutation,
    ):
        self.max_mutation = max_mutation

    def scale(self, depth, a, maxdepth):
        if random() > 0.2:
            pass
        else:
            a = 0

        depthscale = np.exp(-a * (depth / maxdepth) ** 2)
        return depthscale

    def perturb(self, structData, depth, a, maxdepth):

        x = structData["parameters"].copy()
        species = structData["species"].copy()
        u = np.random.normal(0.0, 1.0, len(x))
        delta = (
            (u / np.linalg.norm(u)) * self.max_mutation * self.scale(depth, a, maxdepth)
        )
        x += delta
        x[x > 1] = 1
        x[x < 0] = 0

        return {"parameters": x, "species": species}

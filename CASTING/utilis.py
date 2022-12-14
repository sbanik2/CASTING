"""
Created on 2022-12-13 20:36:01.409428
@author: suvobanik
"""


import numpy as np
import pandas as pd
from pymatgen.core import Lattice


def r_datafame(distance_dict):
    species = []
    for key in distance_dict.keys():
        col = key.split("-")
        species += col
    species = list(set(species))
    data = np.empty((len(species), len(species)))

    for i, sp1 in enumerate(species):
        for j, sp2 in enumerate(species):

            try:
                data[i, j] = distance_dict[sp1 + "-" + sp2]
            except:
                data[i, j] = distance_dict[sp2 + "-" + sp1]

    df = pd.DataFrame(data, index=species, columns=species)

    return df


def get_lattice(box_dim):
    return Lattice.from_parameters(
        a=box_dim, b=box_dim, c=box_dim, alpha=90, beta=90, gamma=90
    )


def get_factors(n, m):

    factors = []
    factor = int(
        n ** (1.0 / m) + 0.1
    )  # fudged to deal with precision problem with float roots
    while n % factor != 0:
        factor = factor - 1
    factors.append(factor)
    if m > 1:
        factors = factors + get_factors(n / factor, m - 1)

    return factors


def DistanceMatrix(frac_coordinates, M):

    a, b, c = (
        frac_coordinates[:, 0],
        frac_coordinates[:, 1],
        frac_coordinates[:, 2],
    )

    def getDist(mat):
        n, m = np.meshgrid(mat, mat)
        dist = m - n
        dist -= np.rint(dist)
        return dist

    da, db, dc = (
        getDist(a),
        getDist(b),
        getDist(c),
    )  # Fracrtional difference matrix

    # ---------cartesian differences------------

    DX = M[0][0] * da + M[1][0] * db + M[2][0] * dc
    DY = M[0][1] * da + M[1][1] * db + M[2][1] * dc
    DZ = M[0][2] * da + M[1][2] * db + M[2][2] * dc

    # -----------distance matrix--------------

    D = np.sqrt(np.square(DX) + np.square(DY) + np.square(DZ))

    return D

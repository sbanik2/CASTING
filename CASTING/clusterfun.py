"""
Created on 2022-12-13 20:36:01.409428
@author: suvobanik
"""


from collections import Counter
from itertools import product
from random import choice, random, shuffle

import networkx as nx
import numpy as np
from ase.build import bulk, make_supercell
from pymatgen.core import Structure
from scipy.spatial.distance import pdist, squareform

from CASTING.utilis import DistanceMatrix, get_factors


def random_sub_cluster_sample(A, natoms):
    indexes = np.arange(0, natoms, 1)
    sampled_nodes = [choice(indexes)]
    for i in range(natoms - 1):
        connections = np.where(A[sampled_nodes[i], :] != 0)[0]
        connections_wth_sam = [n for n in connections if n not in sampled_nodes]
        shuffle(connections_wth_sam)
        next_node = connections_wth_sam[0]
        sampled_nodes.append(next_node)

    return sampled_nodes


def createRandomData(constrains, multiplier=10):

    natoms = constrains["atoms"]

    # --------Species list creation -----------

    species = []
    for key in constrains["composition"].keys():
        nkey = int(
            (
                constrains["composition"][key]
                / sum(list(constrains["composition"].values()))
            )
            * natoms
        )
        for _ in range(nkey):
            species.append(key)

    shuffle(species)

    # ----------------Bulk FCC with natoms X multiplier atoms and select clsuter with random walk -----------

    maxmin = constrains["r_min"].values.max()
    minmax = constrains["r_max"].values.min()
    r = (
        maxmin + random() * abs(minmax - maxmin)
    ) * 0.5  # radius of an indivisual atom from fcc packing
    a = 2 * 2**0.5 * r
    a1 = bulk("Cu", "fcc", a=a)  # 'Cu' dummy species
    factors = get_factors(natoms * multiplier, 3)
    P = np.diag(factors, k=0)  # transformation mattrx
    pos = make_supercell(a1, P).get_positions()
    A = squareform(pdist(pos, metric="euclidean"))
    A[A > 2 * r] = 0
    A[A != 0] = 1
    while True:
        try:
            sampled_nodes = random_sub_cluster_sample(A, natoms)
            break
        except:
            continue
    M = constrains["lattice"].matrix  # lattice matrix
    cluster_pos = pos[sampled_nodes, :]
    box_centre = np.sum(M, axis=0) * 0.5
    cluster_centre = np.mean(cluster_pos, axis=0)
    cluster_pos = box_centre + cluster_centre - cluster_pos
    cluster_pos_fractional = np.matmul(cluster_pos, np.linalg.inv(M))

    return {"parameters": cluster_pos_fractional.flatten(), "species": species}


def get_coords(parameters):
    coords = parameters
    coords = coords.reshape(int(coords.shape[0] / 3), 3)
    return coords


def check_constrains(structData, constrains, verbose=False):

    parameters = structData["parameters"].copy()
    species = structData["species"].copy()
    specieCount = dict(Counter(species))
    coords = get_coords(parameters)

    # ======= number of atom constrains==============

    natoms = constrains["atoms"]
    if natoms != coords.shape[0]:
        if verbose:
            print("# of atoms inconsistent.")
        return False

    # ============ composition check=================

    composition = constrains["composition"]

    if list(composition.keys()).sort() != list(specieCount.keys()).sort():
        if verbose:
            print("composition inconsistent.")
        return False

    for key in composition.keys():
        if specieCount[key] % composition[key] != 0:
            if verbose:
                print("composition inconsistent.")
            return False

    # ======================================

    latt = constrains["lattice"]
    M = np.array((latt.matrix))
    D = DistanceMatrix(coords, M)
    np.fill_diagonal(D, 1e300)
    A = np.empty(D.shape)

    indices = {}

    for key in constrains["composition"].keys():
        indices[key] = np.where(np.array(species) == key)[0].tolist()

    for c in product(list(indices.keys()), repeat=2):

        index1, index2 = indices[c[0]], indices[c[1]]
        d_sub = D[index1, :][:, index2]

        if (d_sub < constrains["r_min"].loc[c[0], c[1]]).any():  # overlapping test
            if verbose:
                print("overlapping atoms.")
            return False

        d_sub[d_sub <= constrains["r_max"].loc[c[0], c[1]]] = 1
        d_sub[d_sub > constrains["r_max"].loc[c[0], c[1]]] = 0

        C = A[index1, :]
        C[:, index2] = d_sub
        A[index1, :] = C

    G = nx.from_numpy_matrix(A)

    # -------fragmentation test----------------

    if len(list(nx.connected_components(G))) > 1:
        if verbose:
            print("Fragmented cluster.")
        return False

    return True


# ---------------------------------------


def parm2struc(structData, constrains):

    parameters = structData["parameters"].copy()
    species = structData["species"].copy()
    pos = get_coords(parameters)
    lattice = constrains["lattice"]
    struct = Structure(lattice, species, pos, to_unit_cell=True)

    return struct


# ----------------------------------------


def struc2param(struct, energy, constrains, CheckFrConstrains=False, writefile=None):

    lattice = constrains["lattice"]
    pos = np.array([list(site.frac_coords) for site in struct.sites]).flatten()

    species = [site.specie.symbol for site in struct.sites]
    latt = [
        lattice.a,
        lattice.b,
        lattice.c,
        lattice.alpha,
        lattice.beta,
        lattice.gamma,
    ]

    StructData = {"parameters": pos, "species": species}

    if CheckFrConstrains:
        if check_constrains(StructData, constrains, verbose=False):
            pass
        else:
            return StructData, 1e300

    if writefile is not None:

        dataString = (
            " ".join(map(str, latt + pos.tolist()))
            + "|"
            + " ".join(species)
            + "|"
            + "{}".format(energy)
        )

        with open(writefile, "a") as outfile:
            outfile.write("{}\n".format(dataString))

    return StructData, energy

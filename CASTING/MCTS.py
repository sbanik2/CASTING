"""
Created on 2022-12-13 20:36:01.409428
@author: suvobanik
"""


from math import log

import numpy as np


def playouts(
    idx,
    depthlist,
    playoutdata,
    Scorelist,
    parameterlist,
    perturbate,
    evaluate,
    a,
    maxdepth,
    nplayouts=10,
):

    nodeID = idx

    for i in range(nplayouts):
        idx += 1
        playdata = perturbate(
            parameterlist[nodeID], depth=depthlist[nodeID], a=a, maxdepth=maxdepth
        )
        playdata_relaxed, playscore = evaluate(playdata)
        #        print("Node: {}, Playout: {} Score: {}".format(nodeID,i+1,playscore))

        Scorelist.append(playscore)
        parameterlist.append(playdata_relaxed)
        playoutdata[nodeID].append(idx)

    return idx, playoutdata, Scorelist, parameterlist


def expansion_simulation(
    parentID,
    idx,
    visits,
    childlist,
    parent,
    depthlist,
    playoutdata,
    Scorelist,
    parameterlist,
    perturbate,
    evaluate,
    a,
    maxdepth,
    nplayouts=10,
):

    playindexes = playoutdata[parentID]
    playscores = [Scorelist[i] for i in playindexes]
    bestplayindex = playindexes[playscores.index(min(playscores))]

    # =========== update =============

    idx += 1
    visits[parentID] += 1
    visits[idx] = 1
    try:
        childlist[parentID].append(idx)
    except:
        childlist[parentID] = [idx]

    parent[idx] = parentID
    depth = depthlist[parentID] + 1
    depthlist[idx] = depth
    playoutdata[idx] = []
    data = perturbate(parameterlist[bestplayindex], depth=depth, a=a, maxdepth=maxdepth)
    data_relaxed, score = evaluate(data)
    Scorelist.append(score)
    parameterlist.append(data)

    idx, playoutdata, Scorelist, parameterlist = playouts(
        idx,
        depthlist,
        playoutdata,
        Scorelist,
        parameterlist,
        perturbate,
        evaluate,
        a,
        maxdepth,
        nplayouts=nplayouts,
    )

    return (
        visits,
        childlist,
        parent,
        depthlist,
        playoutdata,
        Scorelist,
        parameterlist,
        idx,
    )


def backpropagation_selection(
    visits,
    childlist,
    parent,
    depthlist,
    playoutdata,
    Scorelist,
    parameterlist,
    maxdepth,
    exploreconstant=1,
):
    def get_linegae(node, childlist):
        total_list = []
        try:
            total_list += childlist[node]
        except:
            return total_list
        for child in childlist[node]:
            total_list += get_linegae(child, childlist)

        return total_list

    selectionScores = []
    nodes = list(parent.keys())

    for n in nodes:
        par = parent[n]
        depth = depthlist[n]
        if par is None or depth > maxdepth:
            selectionScores.append(-1e300)

        else:
            lineage = get_linegae(n, childlist)
            allindexes = []
            for l in [n] + lineage:
                allindexes += playoutdata[l]

            bestreward = min([Scorelist[i] for i in lineage + allindexes])

            UCB_score = -bestreward + exploreconstant * np.sqrt(
                log(visits[par]) / visits[n]
            )

            selectionScores.append(UCB_score)

    c = zip(selectionScores, nodes)
    c = sorted(c, key=lambda x: -x[0])
    selectionScores, nodes = zip(*c)

    return nodes[0]


def MCTS(
    rootdata,
    perturbate,
    evaluate,
    niterations=200,
    headexpand=10,
    nexpand=3,
    nsimulate=3,
    nplayouts=10,
    exploreconstant=1,
    maxdepth=12,
    a=3,
    selected_node=0,
):

    if selected_node == 0:
        visits = {0: 1}
        childlist = {0: []}
        parent = {0: None}
        depthlist = {0: 0}
        playoutdata = {0: []}
        data_relaxed, score = evaluate(rootdata)
        Scorelist = [score]
        parameterlist = [data_relaxed]
        idx = 0

        # =======run playouts for rootnode ==================

        idx, playoutdata, Scorelist, parameterlist = playouts(
            idx,
            depthlist,
            playoutdata,
            Scorelist,
            parameterlist,
            perturbate,
            evaluate,
            a,
            maxdepth,
            nplayouts=nplayouts,
        )

    # =======simulation and expansion==================

    for _ in range(niterations):

        for _ in range(nexpand):
            (
                visits,
                childlist,
                parent,
                depthlist,
                playoutdata,
                Scorelist,
                parameterlist,
                idx,
            ) = expansion_simulation(
                selected_node,
                idx,
                visits,
                childlist,
                parent,
                depthlist,
                playoutdata,
                Scorelist,
                parameterlist,
                perturbate,
                evaluate,
                a,
                maxdepth,
                nplayouts=nplayouts,
            )

            selected_node = backpropagation_selection(
                visits,
                childlist,
                parent,
                depthlist,
                playoutdata,
                Scorelist,
                parameterlist,
                maxdepth,
                exploreconstant=exploreconstant,
            )

        print(
            "total evaluations: {}, best score yet : {}".format(
                len(Scorelist), min(Scorelist)
            )
        )

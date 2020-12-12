import numpy as np
from functools import partial
from circuitSimHelpers import *
from circuitSimParts import *
from scipy.optimize import newton
import time
import copy

t0 = time.process_time()
##### Options #####
spiceRefNode = '0'
capAdmittance = 0
solverTolerance = 0.0001
refV = 0.0
printRead = True
printResults = True
printSupernodes = True

##### File Input #####
elemDict, nodeDict, _ = circuitPreprocess('./circuit.txt')
if printRead:
    print('**** Preprocessing ****')
    print('Elements:')
    for eid in elemDict.keys():
        elemDict[eid].prnt()
    print('Nodes:')
    print(nodeDict)
else:
    print('**** File Input Complete - Node and Element Definitions Created ****')

##### Net Pre-Processing #####
print('**** Pre-Processing of Net ****')
nodeDictL = copy.deepcopy(nodeDict)
# open circuits:  delete element references from nodeDict
elemstoOpen = set()
for elemid in elemDict.keys():
    if elemDict[elemid].typ == 'C':
        elemstoOpen.add(elemid)
        elemstoOpen.add(-elemid)
for nodeid in nodeDictL.keys():
    nodeDictL[nodeid].elemSet = nodeDictL[nodeid].elemSet - elemstoOpen
# big nodes:
nodeClumps = []
nodesInsideSupernodes = set()
for elemid in elemDict.keys():
    elem = elemDict[elemid]
    if elem.pnode == elem.nnode:
        print('shorted elements not dealt with yet')
        nodeClumps.append((elem.pnode, elem.nnode))
        nodesInsideSupernodes.add(elem.pnode)
        nodesInsideSupernodes.add(elem.nnode)
    elif elem.typ in {'L', 'V'}:
        nodeClumps.append((elem.pnode, elem.nnode))
        nodesInsideSupernodes.add(elem.pnode)
        nodesInsideSupernodes.add(elem.nnode)
nodeClumps = clumpClumps(nodeClumps)
supernodes = dict()
for i, clump in enumerate(nodeClumps):
    supernodes['S0' + str(i+1)] = Supernode(i+1, clump, nodeDictL, elemDict)
if printSupernodes:
    if len(supernodes) > 0:
        print('Supernodes:')
        for snode in supernodes.values():
            snode.prnt()
    else:
        print('- No Supernodes -')

# each supernode knows:  BV, the node collection, the v(BV) function for each collection node, and all attached elems
# the nodeValues loop needs to change so it updates supernodes but not nodes inside them
# the nodeCurrent function needs to be able to take supernodes as inputs
nodesToUpdate = [nid for nid in sorted(nodeDict.keys()) if nid not in nodesInsideSupernodes]
nodesToUpdate.extend([nid for nid in supernodes.keys()])

# Nodal Voltage Update:
for nodeid in nodeDictL:
    nodeDictL[nodeid].V = refV
# Loop through nodes to update voltages
for i in range(10000):
    vchange = 0
    for nodeid in nodesToUpdate:
        if nodeid[0] == 'S':
            currentCalc = partial(nodeCurrent, elemDict, nodeDictL, supernodes[nodeid].elemSet, nUpdate = supernodes[nodeid].nUpdate)
            supernodes[nodeid].BV = newton(currentCalc, supernodes[nodeid].BV)
            vchange = supernodes[nodeid].updateNodes(nodeDictL)
        else:
            vOld = nodeDictL[nodeid].V
            if nodeid == spiceRefNode:
                pass
            else:
                currentCalc = partial(nodeCurrent, elemDict, nodeDictL, {elem: nodeid for elem in nodeDictL[nodeid].elemSet})
                nodeDictL[nodeid].V = newton(currentCalc, nodeDictL[nodeid].V)
            vchange = max(abs(nodeDictL[nodeid].V - vOld), vchange)

    # Post-Process the net
    for nodeid in nodeDictL.keys():
        nodeDict[nodeid].V = nodeDictL[nodeid].V

    ##### Post-Processing #####
    # Loop through elements to update currents
    for eid in elemDict.keys():
        if elemDict[eid].typ == 'I':
            elemDict[eid].current = elemDict[eid].value
        elif elemDict[eid].typ == 'R':
            elemDict[eid].current = (nodeDict[elemDict[eid].pnode].V - nodeDict[elemDict[eid].nnode].V) / elemDict[eid].value
        elif elemDict[eid].typ == 'C':
            elemDict[eid].current = capAdmittance*(nodeDict[elemDict[eid].pnode].V - nodeDict[elemDict[eid].nnode].V)
    # Elements in supernodes:
    for eid in elemDict.keys():
        if elemDict[eid].typ in {'V', 'L'}:
            elemDict[eid].current = 0
            for attachedElem in nodeDict[elemDict[eid].pnode].elemSet:
                if abs(attachedElem) != elemDict[eid].id:
                    elemDict[eid].current += -np.sign(attachedElem)*elemDict[abs(attachedElem)].current

    if vchange < solverTolerance:
        print('\n**** Solver Finished In ' + str(i) + ' Iterations ****\n')
        break
else:
    print('\n***** Solver Finished Without Converging *****\n')

print('**** Postprocessing and Current Calculations Complete ****')
if printResults:
    print('Node Voltages:')
    for nid in nodeDict.keys():
        print('Node ' + str(nid) + ': ' + str(round(nodeDict[nid].V, 2)))
    print('Element Currents:')
    for eid in elemDict.keys():
        if elemDict[eid].current is None:
            print(elemDict[eid].name + ': None')
        else:
            print(elemDict[eid].name + ': ' + str(round(elemDict[eid].current, 2)))
    print('Time Elapsed: ' + str(round(time.process_time() - t0, 4)) + ' s')

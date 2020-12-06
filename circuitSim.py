import numpy as np
from functools import partial
from circuitSimHelpers import clumpClumps
from circuitSimParts import nodeCurrent
from circuitSimParts import circuitPreprocess
from scipy.optimize import newton
import time
import copy

class Supernode:
    def __init__(self, id, nodeCol, nUpdate, nDict, eDict):
        self.id = id
        self.nodeCol = nodeCol
        self.elemSet = set()
        self.BV = 0
        self.nUpdate = nUpdate
        for nid in self.nodeCol:
                self.elemSet.update(nDict[nid].elemSet)
        # take out any elements that appear twice (ones that are inside the supernode)
        posSet = {item for item in self.elemSet if item > 0}
        negSet = {abs(item) for item in self.elemSet if item < 0}
        elems = self.elemSet - posSet.intersection(negSet) - {-item for item in posSet.intersection(negSet)}
        self.elemSet = dict()
        for eid in elems:
            if eid > 0:
                self.elemSet[eid] = eDict[eid].pnode
            else:
                self.elemSet[eid] = eDict[-eid].nnode

    def updateNodes(self, nDict):
        for node in self.nodeCol:
            nDict[node].V = self.BV + self.nUpdate[node]

    def prnt(self):
        print('Node ' + str(self.id) + ':')
        print('    Nodecol: ' + ''.join([str(n)+', ' for n in self.nodeCol]))
        print('    BV: ' + str(self.BV))
        print('    Elemset: ')
        print(self.elemSet)
        print('    nUpdate:')
        print(self.nUpdate)

t0 = time.process_time()
##### Pre-Processing and Solver Options #####
spiceRefNode = 0
capAdmittance = 0
solverTolerance = 0.0001
refV = 0.0
elemDict, nodeDict, _ = circuitPreprocess('./circuit.txt')

# Output Circuit Info:
print('**** Preprocessing ****')
print('Elements:')
for eid in elemDict.keys():
    elemDict[eid].prnt()
print('Nodes:')
print(nodeDict)

##### Solver #####
# Pre-process the net:  set ground, short inductors, open caps
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
        nodeClumps.append((elem.pnode, elem.nnode))
        nodesInsideSupernodes.add(elem.pnode)
        nodesInsideSupernodes.add(elem.nnode)
    elif elem.typ == 'L':
        nodeClumps.append((elem.pnode, elem.nnode))
        nodesInsideSupernodes.add(elem.pnode)
        nodesInsideSupernodes.add(elem.nnode)
    elif elem.typ == 'V':
        nodeClumps.append((elem.pnode, elem.nnode))
        nodesInsideSupernodes.add(elem.pnode)
        nodesInsideSupernodes.add(elem.nnode)
nodeClumps = clumpClumps(nodeClumps)
supernodes = dict()
for i, clump in enumerate(nodeClumps):
    supernodes['S0' + str(i+1)] = Supernode(i+1, clump, {nid: 0 for nid in clump}, nodeDictL, elemDict)
if len(supernodes) > 0:
    print('Supernodes:')
    for snode in supernodes.values():
        snode.prnt()
else:
    print('- No Supernodes -')
#supernodes['S03'].nUpdate = {'N004': 5, 'N006': 0, 'P004': 0}
#supernodes['S01'].nUpdate = {'N002': 10, 'N011': 0}

# each supernode knows:  BV, the node collection, the v(BV) function for each collection node, and all attached elems
# the nodeValues loop needs to change so it updates supernodes but not nodes inside them
# the nodeCurrent function needs to be able to take supernodes as inputs
nodesToUpdate = [nid for nid in sorted(nodeDict.keys()) if nid not in nodesInsideSupernodes]
nodesToUpdate.extend([nid for nid in supernodes.keys()])

# Nodal Voltage Update:
for nodeid in nodeDictL:
    nodeDictL[nodeid].V = refV
# Loop through nodes to update voltages
for i in range(100):
    vchange = 0
    for nodeid in nodesToUpdate:
        if nodeid[0] ==  'S':
            currentCalc = partial(nodeCurrent, elemDict, nodeDictL, supernodes[nodeid].elemSet, nUpdate = supernodes[nodeid].nUpdate)
            supernodes[nodeid].BV = newton(currentCalc, supernodes[nodeid].BV)
            supernodes[nodeid].updateNodes(nodeDictL)
        else:
            vOld = nodeDictL[nodeid].V
            if nodeid == spiceRefNode:
                pass
            else:
                currentCalc = partial(nodeCurrent, elemDict, nodeDictL, {elem: nodeid for elem in nodeDictL[nodeid].elemSet})
                nodeDictL[nodeid].V = newton(currentCalc, nodeDictL[nodeid].V)
            vchange = max(abs(nodeDictL[nodeid].V - vOld), vchange)
    if vchange < solverTolerance:
        print('\n**** Solver Finished In ' + str(i) + ' Iterations ****\n')
        break
else:
    print('\n***** Solver Finished Without Converging *****\n')

# Post-Process the net
for nodeid in nodeDictL.keys():
    nodeDict[nodeid].V = nodeDictL[nodeid].V

##### Post-Processing #####
# Loop through elements to update currents
doCurrents = True
if doCurrents:
    for eid in elemDict.keys():
        if elemDict[eid].typ == 'I':
            elemDict[eid].current = elemDict[eid].value
        elif elemDict[eid].typ == 'R':
            elemDict[eid].current = (nodeDict[elemDict[eid].pnode].V - nodeDict[elemDict[eid].nnode].V) / elemDict[eid].value
        elif elemDict[eid].typ == 'C':
            elemDict[eid].current = capAdmittance*(nodeDict[elemDict[eid].pnode].V - nodeDict[elemDict[eid].nnode].V)
    for eid in elemDict.keys():
        if elemDict[eid].typ == 'V' or elemDict[eid].typ == 'L':
            elemDict[eid].current = 0
            for attachedElem in nodeDict[elemDict[eid].pnode].elemSet:
                if abs(attachedElem) != elemDict[eid].id:
                    elemDict[eid].current += -np.sign(attachedElem)*elemDict[abs(attachedElem)].current

print('**** Postprocessing ****')
print('Node Voltages:')
for nid in nodeDict.keys():
    print('Node ' + str(nid) + ': ' + str(round(nodeDict[nid].V, 2)))
if doCurrents:
    print('Element Currents:')
    for eid in elemDict.keys():
        if elemDict[eid].current is None:
            print(elemDict[eid].name + ': None')
        else:
            print(elemDict[eid].name + ': ' + str(round(elemDict[eid].current, 2)))
print('Time Elapsed: ' + str(round(time.process_time() - t0, 4)) + ' s')

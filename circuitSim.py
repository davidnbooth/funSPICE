import numpy as np
from functools import partial
import sys
sys.path.insert(0, '..')
import funSPICE.circuitSimHelpers as helper
import funSPICE.circuitSimParts as parts
from scipy.optimize import newton
import time
import copy


def funSPICE(inputFile, solverOptions, outputOptions):
    #print('\n\n\n\n**** Simulation Start: ' + inputFile + ' ****')
    # Debugging:
    nodestocheck = set()  # {'N009'}

    spiceRefNode = solverOptions['spiceRefNode']
    capAdmittance = solverOptions['capAdmittance']
    solverTolerance = solverOptions['solverTolerance']
    refV = solverOptions['refV']
    wrelax = solverOptions['wrelax']
    printRead = outputOptions['printRead']
    printResults = outputOptions['printResults']
    printSupernodes = outputOptions['printSupernodes']
    ##### File Input #####
    elemDict, nodeDict, shortedElems, _ = parts.circuitPreprocess(inputFile)
    if printRead:
        print('**** Preprocessing ****')
        print('Elements:')
        for eid in elemDict.keys():
            elemDict[eid].prnt()
        print('Nodes:')
        print(nodeDict)
        for ename in {e for e in shortedElems.values()}:
            print('Element ' + ename + ' is shorted and will be ignored')
    else:
        print('**** File Input Complete - Node and Element Definitions Created ****')

    t0 = time.process_time()
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
        if elem.typ in {'L', 'V'}:
            nodeClumps.append((elem.pnode, elem.nnode))
            nodesInsideSupernodes.add(elem.pnode)
            nodesInsideSupernodes.add(elem.nnode)
    nodeClumps = helper.clumpClumps(nodeClumps)
    supernodes = dict()
    groundClump = None
    for clump in nodeClumps:
        if '0' in clump:
            groundClump = clump
            nodeClumps = [item for item in nodeClumps if item is not clump]
            break
    if groundClump is not None:
        supernodes['S00'] = parts.Supernode(0, groundClump, nodeDictL, elemDict)
    for i, clump in enumerate(nodeClumps):
        supernodes['S0' + str(i+1)] = parts.Supernode(i+1, clump, nodeDictL, elemDict)
    if printSupernodes:
        if len(supernodes) > 0:
            print('Supernodes:')
            for snode in supernodes.values():
                snode.prnt()
        else:
            print('- No Supernodes -')

    # each supernode knows: BV, the node collection, the v(BV) function for each collection node, and all attached elems
    # the nodeValues loop needs to change so it updates supernodes but not nodes inside them
    # the nodeCurrent function needs to be able to take supernodes as inputs
    nodesToUpdate = [nid for nid in sorted(nodeDict.keys()) if nid not in nodesInsideSupernodes]
    nodesToUpdate.extend([nid for nid in supernodes.keys()])

    # Nodal Voltage Update:
    for nodeid in nodeDictL:
        nodeDictL[nodeid].V = refV
    # Loop through nodes to update voltages
    for i in range(solverOptions['maxIters']):
        vchange = 0
        for nodeid in nodesToUpdate:
            if nodeid[0] == 'S':
                if nodeid == 'S00':
                    supernodes[nodeid].BV = 0
                    vchange = max(vchange, supernodes[nodeid].updateNodes(nodeDictL))  # updates the supernode voltages
                else:
                    currentCalc = partial(parts.nodeCurrent, elemDict, nodeDictL, supernodes[nodeid].elemSet, nodeid, nUpdate = supernodes[nodeid].nUpdate)
                    supernodes[nodeid].BV = newton(currentCalc, supernodes[nodeid].BV)*wrelax + (1-wrelax)*supernodes[nodeid].BV
                    vchange = max(supernodes[nodeid].updateNodes(nodeDictL), vchange)  # updates the voltages in supernode
            else:
                vOld = nodeDictL[nodeid].V
                if nodeid == spiceRefNode:
                    pass
                else:
                    verbose = False
                    if nodeid in nodestocheck:
                        verbose = True
                    currentCalc = partial(parts.nodeCurrent, elemDict, nodeDictL, {elem: nodeid for elem in nodeDictL[nodeid].elemSet}, nodeid, verbose=verbose)
                    newV = newton(currentCalc, nodeDictL[nodeid].V)
                    nodeDictL[nodeid].V = newV*wrelax + (1-wrelax)*vOld
                vchange = max(abs(nodeDictL[nodeid].V - vOld), vchange)

        # Post-Process the net
        for nodeid in nodeDictL.keys():
            nodeDict[nodeid].V = nodeDictL[nodeid].V

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
                eCP = 0
                eCN = 0
                for attachedElem in nodeDict[elemDict[eid].pnode].elemSet:
                    if abs(attachedElem) != elemDict[eid].id:
                        try:
                            eCP += -np.sign(attachedElem)*elemDict[abs(attachedElem)].current
                        except Exception:
                            raise ValueError('eid: ' + str(eid) + ' failed on current update when adding ' + str(abs(attachedElem)))
                for attachedElem in nodeDict[elemDict[eid].nnode].elemSet:
                    if abs(attachedElem) != elemDict[eid].id:
                        try:
                            eCN += np.sign(attachedElem)*elemDict[abs(attachedElem)].current
                        except Exception:
                            raise ValueError('eid: ' + str(eid) + ' failed on current update when adding ' + str(abs(attachedElem)))
                elemDict[eid].current = 0.5*(eCP + eCN)

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
        elems = helper.customListSort(list(elemDict.keys()), (20, 51, 30, 40, 10, 50))
        for eid in elems:
            if elemDict[eid].current is None:
                print(elemDict[eid].name + ': None')
            else:
                print(elemDict[eid].name + ': ' + str(round(elemDict[eid].current, 2)))
        for ename in {e for e in shortedElems.values()}:
            print('Element ' + ename + ' has no current because it was shorted and ignored')
        print('Solver Time Elapsed: ' + str(round(time.process_time() - t0, 4)) + ' s')
    return nodeDict, elemDict, shortedElems


if __name__ == '__main__':
    inputFile = './testcircuits/diodemadness.txt'
    # inputFile = './circuit.txt'
    solverOptions = dict(spiceRefNode='0', capAdmittance=0, solverTolerance=10e-6, refV=0.0, maxIters=int(3*10e4), wrelax=1)
    outputOptions = dict(printRead=True, printResults=True, printSupernodes=True)
    nodeDict, elemDict, shortedElems = funSPICE(inputFile, solverOptions, outputOptions)
    parts.writeOutput(nodeDict, elemDict, './output.txt')

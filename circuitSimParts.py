import numpy as np

def circuitPreprocess(filepath):
    class Elem:
        def __init__(self, name, pnode, nnode, value):
            try:
                self.typ = name[0].upper()
                self.id = {'R': 10000, 'C': 20000, 'L': 30000, 'I': 40000, 'V': 50000}[self.typ] + int(name[1:])
            except Exception:
                print('Element Identifier not understood: ' + self.name)
            # strip leading Ns from node numbers
            self.nnode = nnode
            self.pnode = pnode
            self.value = float(value)
            self.name = name
            self.current = None  # for postprocessing

        def prnt(self):
            print(str(self.name) + ' ' + str(self.pnode) + ' ' + str(self.nnode) + ' ' + str(self.value))

    class Node:
        def __init__(self, id):
            self.id = id
            self.elemSet = set()
            self.V = 0

    # Reading in the file with comments
    with open(filepath, 'r') as f:
        inputCirc = []
        for line in f:
            if line[0][0] not in {'.', '%', '*', '#'}:
                inputCirc.append(line.split())
    inputCirc = np.array(inputCirc, dtype=str)

    # Make Node and Element Dicts
    elemDict = dict()  # element ids are keys, value is element object
    nodeDict = dict()  # node ids are keys, value is node object.  in the node object eIDs are (-) if node is at (-) end
    vsources = []  # pnode, nnode, value, typ, id for sources
    for row in inputCirc:
        if row[0][0] not in {'.', '%', '#'}:
            # add to element dictionary:
            elem = Elem(*row)
            elemDict[elem.id] = elem
            # add to node dictionary:
            pnode = elem.pnode
            nnode = elem.nnode
            id = elem.id
            if pnode in nodeDict.keys():
                nodeDict[pnode].elemSet.add(id)
            else:
                nodeDict[pnode] = Node(pnode)
                nodeDict[pnode].elemSet.add(id)
            if nnode in nodeDict.keys():
                nodeDict[nnode].elemSet.add(-id)
            else:
                nodeDict[nnode] = Node(nnode)
                nodeDict[nnode].elemSet.add(-id)
            if elem.typ == 'V':
                vsources.append([pnode, nnode, float(elem.value), elem.typ, elem.id])

    return (elemDict, nodeDict, vsources)

def nodeCurrent(elemDict, nodeDictL, attachedElems, startV, nUpdate=None):
    nodalCurrent = 0
    for eID, nid in attachedElems.items():
        if nUpdate is None:
            V = startV
        else:
            V = startV + nUpdate[nid]
        if elemDict[abs(eID)].typ == 'R':
            if eID > 0:
                resistorV = nodeDictL[elemDict[eID].nnode].V - V
            else:
                resistorV = nodeDictL[elemDict[abs(eID)].pnode].V - V
            nodalCurrent += resistorV / elemDict[abs(eID)].value
        elif elemDict[abs(eID)].typ == 'I':
            nodalCurrent += -elemDict[abs(eID)].value * np.sign(eID)
        else:
            print('son u have a problem in the nodeCurrent calc')
            print(elemDict[abs(eID)].typ)
    return nodalCurrent

import numpy as np
from circuitSimHelpers import supernodeInternalGraph


class Elem:
    def __init__(self, name, pnode, nnode, value):
        eTypeDict = {'R': 10000, 'C': 20000, 'L': 30000, 'I': 40000, 'V': 50000}
        try:
            self.typ = name[0].upper()
            self.id = eTypeDict[self.typ] + int(name[1:])
            self.name = name
        except Exception:
            print('Element Identifier not understood: ' + self.name)
        self.nnode = nnode
        self.pnode = pnode
        self.value = float(value)
        self.current = None

    def prnt(self):
        print(str(self.name) + ' ' + str(self.pnode) + ' ' + str(self.nnode) + ' ' + str(self.value))


class Node:
    def __init__(self, id):
        self.id = id
        self.elemSet = set()
        self.V = 0


class Supernode:
    def __init__(self, id, nodeCol, nDict, eDict):
        self.id = id
        self.nodeCol = nodeCol
        self.elemSet = set()
        self.BV = 0
        for nid in self.nodeCol:
                self.elemSet.update(nDict[nid].elemSet)
        # find elements internal to the supernode:
        self.internalElemset = set()
        for eid in self.elemSet:
            if -eid in self.elemSet and eDict[abs(eid)].typ in {'V', 'L'}:
                self.internalElemset.add(abs(eid))

        elems = self.elemSet - self.internalElemset - {-item for item in self.internalElemset}
        self.elemSet = dict()
        for eid in elems:
            if eid > 0:
                self.elemSet[eid] = eDict[eid].pnode
            else:
                self.elemSet[eid] = eDict[-eid].nnode
        #define internalNodes - a list of node objects without references to elems outside the supernode and also internalElems
        self.internalNodes = dict()
        for nid in self.nodeCol:
            self.internalNodes[nid] = Node(nid)
            self.internalNodes[nid].elemSet = nDict[nid].elemSet.intersection(self.internalElemset) | nDict[nid].elemSet.intersection({-item for item in self.internalElemset})
            self.internalNodes[nid].V = None
        self.internalElems = dict()
        for eid in self.internalElemset:
            elem = Elem(eDict[eid].name, eDict[eid].pnode, eDict[eid].nnode, 0)
            self.internalElems[eid] = elem
            if elem.typ == 'V':
                self.internalElems[eid].value = eDict[eid].value
        self.nUpdate = supernodeInternalGraph(self.internalNodes, self.internalElems)


    def updateNodes(self, nDict):
        vchange = 0
        for node in self.nodeCol:
            oldV = nDict[node].V
            nDict[node].V = self.BV + self.nUpdate[node]
            vchange = max(vchange, abs(oldV - nDict[node].V))
        return vchange

    def prnt(self):
        print('Node ' + str(self.id) + ':')
        print('    Nodecol: ' + ''.join([str(n)+', ' for n in self.nodeCol]))
        print('    BV: ' + str(self.BV))
        print('    Elemset: ')
        print('            ' + ''.join([str(k) + ': ' + str(v) + ', ' for k, v in self.elemSet.items()]))
        print('    nUpdate:')
        print('            ' + ''.join([str(k) + ': ' + str(v) + ', ' for k, v in self.nUpdate.items()]))
        print('    Internal Elements:')
        for elem in self.internalElems.values():
            print('        ' + str(elem.id) + ':  ' + str(elem.pnode) + ' ' + str(elem.nnode))
        print('    Internal Nodes:')
        for node in self.internalNodes.values():
            print('        ' + str(node.id) + ':  ' + ''.join([str(eid)+', ' for eid in node.elemSet]))


def circuitPreprocess(filepath):
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

    # Take out any shorted elements
    shortedElems = set()
    for eid in elemDict.keys():
        if elemDict[eid].pnode == elemDict[eid].nnode:
            shortedElems.add(abs(eid))
            shortedElems.add(-abs(eid))
    for eid in shortedElems:
        if eid > 0:
            del elemDict[eid]
    for node in nodeDict.values():
        node.elemSet = node.elemSet - shortedElems

    return (elemDict, nodeDict, shortedElems, vsources)

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

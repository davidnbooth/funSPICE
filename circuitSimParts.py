import numpy as np
from circuitSimHelpers import supernodeInternalGraph
import copy


class Elem:
    def __init__(self, name, pnode, nnode, value):
        eTypeDict = {'R': 10000, 'C': 20000, 'L': 30000, 'I': 40000, 'V': 50000}
        self.equation = dict()
        try:
            if name[0].upper() == 'D':
                self.typ = 'V'
                self.name = name
                self.equation['n'] = 1
                self.equation['Is'] = 10e-14
                self.id = eTypeDict[self.typ] + int(name[1:]) + 1000
            else:
                self.typ = name[0].upper()
                self.id = eTypeDict[self.typ] + int(name[1:])
                self.name = name
        except Exception:
            print('Element Identifier not understood: ' + self.name)
        self.nnode = nnode
        self.pnode = pnode
        if len(self.equation) > 0:
            self.value = -0.8
        else:
            self.value = float(value)
        if self.typ in {'V', 'L'}:
            self.current = 1
        else:
            self.current = None

    def updateValue(self):
        # Diode:
        if len(self.equation) > 0:
            n = self.equation['n']
            Vt = 0.02585
            Is = self.equation['Is']
            self.value = n*Vt*np.log(self.current/Is+1)

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
        # define internalNodes - a list of node objects without references to elems outside the supernode and also internalElems
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
        if self.id == 0:
            refV = self.nUpdate['0']
            for key in self.nUpdate.keys():
                self.nUpdate[key] = self.nUpdate[key] - refV

    def updateNodes(self, nDict):
        for elem in self.internalElems.values():
            elem.updateValue()
        nodesfornupdate = copy.deepcopy(self.internalNodes)
        for node in nodesfornupdate.values():
            node.V = None
        self.nUpdate = supernodeInternalGraph(nodesfornupdate, self.internalElems)
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
    if read_file:
        with open(input, 'r') as f:
            inputCirc = []
            for line in f:
                line = line.strip().replace('\x00', '')
                if line[0][0] not in {'.', '%', '*', '#'}:
                    inputCirc.append(line.split())
        inputCirc = np.array(inputCirc, dtype=str)
    else:
        input = input.strip().split('\n')
        inputCirc = []
        for line in input:
            line = line.strip().replace('\x00', '')
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
    shortedElems = dict()
    for eid in elemDict.keys():
        if elemDict[eid].pnode == elemDict[eid].nnode:
            shortedElems[abs(eid)] = elemDict[abs(eid)].name
            shortedElems[-abs(eid)] = elemDict[abs(eid)].name
    shortedElemIds = {eid for eid in shortedElems.keys()}
    for eid in shortedElemIds:
        if eid > 0:
            del elemDict[eid]
    for node in nodeDict.values():
        node.elemSet = node.elemSet - shortedElemIds

    return elemDict, nodeDict, shortedElems, vsources


def nodeCurrent(elemDict, nodeDictL, attachedElems, nodeid, startV, nUpdate=None, verbose=False):
    nodalCurrent = 0
    # For each attached element calculate its current
    for eID, nid in attachedElems.items():  # attachedElems is {attachedElemID: nodeAttachedtothatelem, ...} -> for supernodes
        if verbose:
            print(elemDict[abs(eID)].name + 'Current: ')
        # define voltage of current node
        if nUpdate is None:
            V = startV
        else:
            V = startV + nUpdate[nid]
        # if its a resistor, calculate the current:
        if elemDict[abs(eID)].typ == 'R':
            if eID > 0:
                resistorV = nodeDictL[elemDict[eID].nnode].V - V
            else:
                resistorV = nodeDictL[elemDict[abs(eID)].pnode].V - V
            resistorI = resistorV / elemDict[abs(eID)].value
            if verbose:
                print(resistorI)
            nodalCurrent += resistorI
        elif elemDict[abs(eID)].typ == 'I':
            sourceI = -elemDict[abs(eID)].value * np.sign(eID)
            nodalCurrent += sourceI
            if verbose:
                print(sourceI)
        else:
            print('son u have a problem in the nodeCurrent calc for node ' + str(nodeid))
            print(elemDict[abs(eID)].typ)
        if verbose:
            print('Voltage at "home" node ' + str(nid) + ': ' + str(round(V, 2)))
    return nodalCurrent


def writeOutput(nDict, eDict, filename):
    with open(filename, 'w') as f:
        for nid in nDict.keys():
            f.write('V(' + str(nid) + ') ' + str(round(nDict[nid].V, 2)) + '\n')
        for eid in eDict.keys():
            if eDict[eid].current is None:
                f.write('I(' + eDict[eid].name + ') None\n')
            else:
                f.write('I(' + eDict[eid].name + ') ' + str(round(eDict[eid].current, 2)) + '\n')
    return True


def readOutput(filename):
    nodeDict = dict()
    elemDict = dict()
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip().replace('\x00', '')
            if line[0] == 'V':
                nodeDict[line.split()[0][2:-1].upper()] = float(line.split()[1])
            elif line[0] == 'I':
                elemDict[line.split()[0][2:-1].upper()] = float(line.split()[1])
    if '0' not in nodeDict.keys():
        nodeDict['0'] = 0.0
    return nodeDict, elemDict


# test file vs ref file
def compareOutputs(tfile, rfile):
    if type(tfile) == str:
        tNodes, tElems = readOutput(tfile)
    else:
        tNodes = tfile[0]
        tElems = tfile[1]
    if type(rfile) == str:
        rNodes, rElems = readOutput(rfile)
    else:
        rNodes = rfile[0]
        rElems = rfile[1]
    nodeDiff = dict()
    for nid in tNodes:
        if nid in rNodes:
            if np.isclose(tNodes[nid], rNodes[nid]):
                nodeDiff[nid] = (0, 0)
            else:
                nodeDiff[nid] = (2 * (tNodes[nid] - rNodes[nid]) / (tNodes[nid] + rNodes[nid]), abs(tNodes[nid] - rNodes[nid]))
        else:
            print('Results are different length!!')
            break
    elemDiff = dict()
    for eid in tElems:
        if eid in rElems:
            if np.isclose(tElems[eid], rElems[eid]):
                elemDiff[eid] = (0, 0)
            else:
                elemDiff[eid] = (2 * (tElems[eid] - rElems[eid]) / (tElems[eid] + rElems[eid]), abs(tElems[eid] - rElems[eid]))
        else:
            print('Results are different length!!')
            break
    return nodeDiff, elemDiff


if __name__ == '__main__':
    print(readOutput('./testcircuits/multiplesourcesRLC_out.txt'))

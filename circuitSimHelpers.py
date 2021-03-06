import numpy as np
#np.random.seed(111)
### BFS for supernodes
#  Each node knows about its edges and has a value
#  Each edge has a directed value and two nodes which it knows about.
#  A set of nodes which have been explored is maintained
#
def supernodeInternalGraph(nodes, edges):
    explored = set()
    queue = []

    # pick a root and define its value as reference (zero)
    root = sorted(list(nodes.keys()))[0]
    queue.append(root)
    nodes[root].V = 0
    # pull a node
    while len(queue) > 0:
        nid = queue.pop()
        explored.add(nid)
        # for each edge, find the attached nodes
        for eid in nodes[nid].elemSet:
            attachedNode = edges[eid].nnode if eid > 0 else edges[abs(eid)].pnode  # get attached node from other side of element
            val = nodes[nid].V - np.sign(eid) * edges[abs(eid)].value  # value attached node should get
            if attachedNode in explored:  # attached node has been fully explored
                if abs(val - nodes[attachedNode].V) > 0.001:
                    print('Invalid Graph')
                    print(attachedNode + ' has value ' + str(nodes[attachedNode].V) + ' and node ' + nid + ' is trying to assign it ' + str(val))
            elif nodes[attachedNode].V is not None:  # attached node has not been visited but has been given a value by another attached node
                if abs(val - nodes[attachedNode].V) > 0.001:
                    print('Invalid Graph')
                    print(attachedNode + ' has value ' + str(nodes[attachedNode].V) + ' and node ' + nid + ' is trying to assign it ' + str(val))
            else:
                queue.append(attachedNode)
                nodes[attachedNode].V = val
    return {nid: nodes[nid].V for nid in nodes.keys()}


def clumpClumps(inputList):     # This code is contributed by Abhishek Valsan
    # Python program to print connected
    # components in an undirected graph
    class Graph:

        # init function to declare class variables
        def __init__(self, V):
            self.V = V
            self.adj = [[] for i in range(V)]

        def DFSUtil(self, temp, v, visited):

            # Mark the current vertex as visited
            visited[v] = True

            # Store the vertex to list
            temp.append(v)

            # Repeat for all vertices adjacent
            # to this vertex v
            for i in self.adj[v]:
                if visited[i] == False:

                    # Update the list
                    temp = self.DFSUtil(temp, i, visited)
            return temp

        # method to add an undirected edge
        def addEdge(self, v, w):
            self.adj[v].append(w)
            self.adj[w].append(v)

        # Method to retrieve connected components
        # in an undirected graph
        def connectedComponents(self):
            visited = []
            cc = []
            for i in range(self.V):
                visited.append(False)
            for v in range(self.V):
                if visited[v] == False:
                    temp = []
                    cc.append(self.DFSUtil(temp, v, visited))
            return cc

    vertices = set()
    for v in inputList:
        vertices.add(v[0])
        vertices.add(v[1])
    g = Graph(len(vertices))
    vertexTransIn = {v: i for i, v in enumerate(vertices)}
    vertexTransOut = dict()
    for k, v in vertexTransIn.items():
        vertexTransOut[v] = k
    for v in inputList:
        g.addEdge(vertexTransIn[v[0]], vertexTransIn[v[1]])
    cc = []
    for sublist in g.connectedComponents():
        cc.append([vertexTransOut[v] for v in sublist])
    return(cc)

#lmao dont need this its in scipy
def secantUpdate(f, initialGuess, maxiterations, tolerance):
    Vn0 = initialGuess
    Vn1 = Vn0 + 1
    for i in range(maxiterations):
        olderfValue = f(Vn0)
        oldfValue = f(Vn1)
        if abs(oldfValue - olderfValue) > tolerance:
            Vn1temp = Vn1
            Vn1 = Vn1 - oldfValue * (Vn1 - Vn0) / (oldfValue - olderfValue)
            Vn0 = Vn1temp
        else:
            break
    return Vn1


# takes a list of equal-digit integers and a list of starts
# e.g.   51001 51000  50000  30001 30000 and (50, 51, 30)
# and returns a list where all numbers that start with 50 are first, then 51, then 30, each sublist is reverse sorted
def customListSort(inList, sortOrder):
    checkDigits = int(np.log10(sortOrder[0]) + 1)
    numDigits = int(np.log10(inList[0]) + 1)
    nums = {n: [] for n in sortOrder}
    for number in inList:
        nums[int(number/(10**(numDigits-checkDigits)))].append(number)
    outputList = []
    for num in sortOrder:
        outputList.extend(reversed(sorted(nums[num])))
    return outputList

if __name__ == '__main__':
    inputList = ((0, 1), (2, 1), (2, 40), (5, 6), (6, 8), (3, 7))
    #print(clumpClumps(inputList))

    #
    #     Example:
    #
    #        +e1-
    #     A ----- B
    #   + |
    #   e2|
    #   - |-e3+   -e4+
    #     C----E----F
    #  +  |   /  -
    #   e6|  / e5
    #  -  | /    +
    #     D
    #
    #  Should be:  A
    exampleNodes = {'A': {1, 2}, 'B': {-1,}, 'C': {-2, -3, 6}, 'D': {5, -6}, 'E': {3, -4, -5}, 'F': {4, }}
    exampleEdges = {'V1': ('A', 'B', 4), 'V2': ('A', 'C', 3), 'V3': ('E', 'C', 2), 'V4': ('F', 'E', 2), 'V5': ('D', 'E', 2), 'V6': ('C', 'D', -4)}
    exampleNodes = {'N002': {50001,}, 'N012': {-50001}}

    from circuitSimParts import Node
    from circuitSimParts import Elem
    nodes = {k: Node(k) for k in exampleNodes.keys()}
    for nid, node in nodes.items():
        node.elemSet = exampleNodes[nid]
    edges = {k: Elem(k, v[0], v[1], v[2]) for k, v in exampleEdges.items()}

    #edges = {50001: Elem('V01', 'N002', 'N012', 10.)}
    #print(supernodeInternalGraph(nodes, edges))

    inList = [51004, 51003, 51002, 51001, 50003, 50002, 50001, 40002, 40001, 30004, 30003, 30002, 30001, 20003, 20001, 10016, 10015, 10014, 10013, 10012, 10011, 10010, 10009, 10008, 10007, 10006, 10005, 10004, 10003, 10002, 10001]
    print(customListSort(inList, (20, 51, 30, 40, 10, 50)))
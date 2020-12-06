
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

if __name__ == '__main__':
    inputList = ((0, 1), (2, 1), (2, 40), (5, 6), (6, 8), (3, 7))
    print(clumpClumps(inputList))
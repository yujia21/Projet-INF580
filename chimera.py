import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

def chimera(k):
    """
    Create a chimera graph(via matrix) for a given k, a chimera graph has 8*k*k vertices
    """
    g = np.zeros((8*k*k,8*k*k))
    # generate a complete bipartite graph K4,4 for k*k grid graph
    for i in range(k*k):
        g[8*i+4:8*i+8,8*i:8*i+4] = np.ones((4,4))
        g[8*i:8*i+4,8*i+4:8*i+8] = np.ones((4,4))
    # right partition
    for i in range(k):
        for j in range(k-1):
            n = 8*k*i+j*8
            g[n+4:n+8,n+12:n+16] = np.diag([1]*4)
            g[n+12:n+16,n+4:n+8] = np.diag([1]*4)
    # left partition
    for i in range(k-1):
        for j in range(k):
            n = 8*k*i+j*8
            g[n:n+4,n+k*8:n+k*8+4] = np.diag([1]*4)
            g[n+k*8:n+k*8+4,n:n+4] = np.diag([1]*4)
    return g

def chimera_nx(k):
    """
    Create a chimera graph(via nx) for a given k
    """
    l=[(0, 4), (0, 5), (0, 6), (0, 7), (1, 4), (1, 5), (1, 6), (1, 7), (2, 4), (2, 5), (2, 6), (2, 7), (3, 4),(3, 5), (3, 6), (3, 7)]
    G = nx.Graph()
    for i in range(8*k*k):
        G.add_node(i)
    # Graph k4,4
    for i in range(k*k):
        G.add_edges_from(map(lambda t:(t[0]+8*i,t[1]+8*i),l))
    for i in range(k):
        for j in range(k):
            n = 8*k*i+8*j
            if not j==k-1:
                G.add_edges_from([(n+4,n+12),(n+5,n+13),(n+6,n+14),(n+7,n+15)])
            if not i==k-1:
                G.add_edges_from([(n,n+8*k),(n+1,n+1+8*k),(n+2,n+2+8*k),(n+3,n+3+8*k)])
    return G

if __name__ == '__main__':
    #g = chimera(2)
    #for i in range(len(g)):
    #    s=''.join(map(str,map(int,g[i])))
    #    print(s)
    g = chimera_nx(2)
    nx.draw(g)
    plt.show()

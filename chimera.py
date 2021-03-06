import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import random

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

def chimera_nx(k,r=0):
    """
    Create a chimera graph(via nx) for a given k, r is number of nodes to remove, optional
    """
    l=[(0, 4), (0, 5), (0, 6), (0, 7), (1, 4), (1, 5), (1, 6), (1, 7), (2, 4), (2, 5), (2, 6), (2, 7), (3, 4),(3, 5), (3, 6), (3, 7)]
    G = nx.Graph()
    pos = {}
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
    for i in range(k):
        for j in range(k):
            p = (15*j,25*i)
            for x in range(4):
                pos[x+8*(j+(i*k))] = (p[0],p[1]+5*x)
            for x in range(4):
                pos[4+x+8*(j+(i*k))] = (p[0]+5,p[1]+5*x)

    remove = random.sample(range(8*k*k),r)
    print(r,remove)
    G.remove_nodes_from(remove)
    G.add_nodes_from(remove)
    return G,pos

def output_chimera_edgelist(k,r) :
    assert r < 8*k*k
    G = chimera_nx(k,r)
    for (i, j) in G.edges() : 
        #random weight for each edge between 0 to 5
        G[i][j]['weight'] = 1 #int(np.ceil(np.random.uniform(0,5)))
    nx.write_weighted_edgelist(G, 'data_chimera.nonweighted.reduced2k2/'+str(k)+'_0.edgelist')

if __name__ == '__main__':
    #g = chimera(2)
    #for i in range(len(g)):
    #    s=''.join(map(str,map(int,g[i])))
    #    print(s)
#    g = chimera_nx(2,0)
#    nx.draw(g)
#    plt.show()
    g,pos = chimera_nx(3,0)
    nx.draw(g,pos,with_labels=True)
    plt.show()
#    for k in [1, 2, 3, 4] :
#        output_chimera_edgelist(k,2*k*k)

import networkx as nx 
import numpy as np

# Creates directed graph with Erdos Renyi Model
# n nodes, p probability for each edge
for n in np.arange(10, 100, 10):
    for p in np.arange(0.3, 0.9, 0.2):
        G=nx.erdos_renyi_graph(n,p)
        for (i, j) in G.edges() : 
            #random weight for each edge between 0 to 5
            G[i][j]['weight'] = int(np.ceil(np.random.uniform(0,5)))
        nx.write_weighted_edgelist(G, 'data_'+str(n)+'_'+str(p)+'.edgelist')
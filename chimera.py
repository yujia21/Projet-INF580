import numpy as np

def chimera(k):
    """
    Create a chimera graph for a given k, a chimera graph has 8*k*k vertices
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

if __name__ == '__main__':
    g = chimera(2)
    for i in range(len(g)):
        s=''.join(map(str,map(int,g[i])))
        print(s)

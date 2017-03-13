import time
import sys
import math
import numpy as np
from pyomo.environ import *
from pyomo.opt import SolverStatus, TerminationCondition
import pyomo.opt
import cvxopt as cvx #used by picos
import picos as pic #can solve SDP
import networkx as nx #manipulate graphs

GWitn = 10000 #max iteration of randomized rounding of GW
MyEpsilon = 1e-6 # tolerance
cplexTimeLimit= 100 #maximum seconds of user CPU for CPLEX

#Quadratic programming : pyomo
def create_maxcut_miqp(G) : 
    maxcut = ConcreteModel()
    # vars = 0 if node in S, 1 if not
    maxcut.y = Var(G.nodes(), within = Binary)
    
    #objective fn
    def max_cut(model):
        return sum(G[i][j]['weight'] * (model.y[i]+model.y[j] - \
                   2*model.y[i]*model.y[j]) for (i,j) in G.edges())
    maxcut.maxcutobj = Objective(rule = max_cut, sense = maximize)
    return maxcut

#Linear programming : pyomo
def create_maxcut_milp(G):
    maxcut = ConcreteModel()
    # vars = 0 if node in S, 1 if not
    maxcut.y = Var(G.nodes(),within = Binary)
    # linearization vars : z(i,j) = y(i)*y(j)
    maxcut.z = Var(G.edges(), within = PercentFraction)
    # objective fn : replace y(i)*y(j) with z(i,j)
    def max_cut(model):
        return sum(G[i][j]['weight'] * (model.y[i]+model.y[j] - \
                   2*model.y[i]*model.z[i,j]) for (i,j) in G.edges())
    maxcut.maxcutobj = Objective(rule = max_cut, sense = maximize)
    # Fortet's linearization constraints : 
    def fortet1(model, i, j):
        return model.z[i,j] <= model.y[i]
    maxcut.fortet1 = Constraint(G.edges(), rule = fortet1)
    def fortet2(model, i, j):
        return model.z[i,j] <= model.y[j]
    maxcut.fortet2 = Constraint(G.edges(), rule = fortet2)
    def fortet3(model, i, j):
        return model.z[i,j] >= model.y[i] + model.y[j] -1
    maxcut.fortet3 = Constraint(G.edges(), rule = fortet3)
    return maxcut
    
#SDP relaxation : picos
def create_maxcut_sdp(Laplacian):
    n = Laplacian.shape[0]
    maxcut = pic.Problem()
    X = maxcut.add_variable('X', (n,n), 'symmetric')
    maxcut.add_constraint(pic.tools.diag_vect(X) == 1)
    maxcut.add_constraint(X >> 0)
    L = pic.new_param('L', Laplacian)
    # max trace(LX)
    maxcut.set_objective('max', L|X)
    return maxcut

#Matrix factor V of X*    
def rank_factor(X, rk):
    n = X.shape[0]
    evl, evc = np.linalg.eigh(X)
    # set neg eigv to 0
    evl[evl<0] = 0
    if rk < n:
        V = np.transpose(evc[:,-rk:])
        for k in range(rk):
            V[k] *= math.sqrt(evl[-rk-1])
    else : 
        V = np.transpose(evc[::-1])
        for k in range(n):
            V[k] *= math.sqrt(evl[::-1][k])
    return V

# Gets closest sdp from X (force eigen > 0)    
def closest_sdp(X) : 
    evl, evc = np.linalg.eigh(X)
    evl[evl < 0] = 0
    return evc.T.dot(np.diag(evl).dot(evc))
    
def gw_randomized_rounding(V, objX, L, iterations):
    n = V.shape[0]
    count = 0
    obj = 0
    print("maxcut(GW/rnd) : impro_ithn objfun")
    while (count<iterations):
        # normalized vector on unit sphere
        r = np.random.normal(0,1,n)
        r /= np.linalg.norm(r) 
        
        # sign of V.r
        x = np.sign(np.dot(V,r))
        x[x>=0] = 1 # avoid 0
        
        # obj fn
        o = np.dot(x.T, np.dot(L,x))
        
        if o > obj : 
            x_cut = x
            obj = o
            print(str(count)+"\t"+str(obj))
        count +=1
    return (x_cut, obj)
    

# MAIN
# read edge list file as input
if len(sys.argv) < 3:
    print("error : need an edge list filename on cmd line and specify if chimera")
    quit()
# each line : i j w[i,j]
G = nx.read_weighted_edgelist(sys.argv[1], nodetype = int)
G = nx.Graph(G)

# MIQP formulation with Pyomo
miqp = create_maxcut_miqp(G)
solver = pyomo.opt.SolverFactory('cplexamp', solver_io = 'nl')
solver.options['timelimit'] = cplexTimeLimit
solver.options['mipdisplay'] = 2
t0 = time.time()
results= solver.solve(miqp, keepfiles = False, tee = True)
t1 = time.time()
miqpcpu = t1 - t0
miqp.solutions.load_from(results)
ymiqp = {i:miqp.y[i].value for i in G.nodes()}
miqpcut = [i for i, y in ymiqp.iteritems() if y == 1]
miqpobj = miqp.maxcutobj()

# MILP formulation
milp = create_maxcut_milp(G)
t0 = time.time()
results = solver.solve(milp, keepfiles = False, tee = True)
t1 = time.time()
milpcpu = t1 - t0
milp.solutions.load_from(results)
ymilp = {i:milp.y[i].value for i in G.nodes()}
milpcut = [i for i, y in ymilp.iteritems() if y == 1]
milpobj = milp.maxcutobj()

# SDP relaxation
L = np.array(1/4.*nx.laplacian_matrix(G).todense())
sdp = create_maxcut_sdp(L)
t0 = time.time()
sdp.solve(verbose = 0)
t1 = time.time()
sdpcpu = t1 - t0
sdpobj = sdp.obj_value()
Xsdp = closest_sdp(np.array(sdp.get_valued_variable('X')))
sdprank = np.linalg.matrix_rank(Xsdp)
sdpfull = Xsdp.shape[0]

#t0 = time.time() # because includes SDP time
N = len(G.nodes())
V = rank_factor(Xsdp, N)
(gwcut, gwobj) = gw_randomized_rounding(V, sdp.obj_value(), L, GWitn)
t1 = time.time()
gwcpu = t1 - t0

if (sys.argv[2] == "erdos_renyi") : 
    # Erdos Renyi Graphs
    n_p_list = sys.argv[1].split('_')
    n = n_p_list[0][-2:]
    p = n_p_list[1][:3]

    print(str(n)+" nodes, with proba "+str(p)+' of any edge')
    print("maxcut(out) : MIQPobj = "+str(miqpobj)+", MIQPcpu ="+str(miqpcpu))
    print("maxcut(out) : MILPobj = "+str(milpobj)+", MILPcpu ="+str(milpcpu))
    print("maxcut(out) : SDPobj = "+str(sdpobj)+", SDPcpu ="+str(sdpcpu))
    print("maxcut(out) : GWobj = "+str(gwobj)+", GWcpu ="+str(gwcpu))
    print('')

    f = open('results.out', 'a')
    f.write(str(n)+' '+str(p)+' ')
    f.write(str(miqpobj)+' '+str(milpobj)+' '+str(sdpobj)+'('+str(sdprank)+'/'+str(sdpfull)+') '+str(gwobj)+' ')
    f.write(str(miqpcpu)+' '+str(milpcpu)+' '+str(sdpcpu)+' '+str(gwcpu)+'\n')
    f.close()

elif (sys.argv[2] == "chimera") :
    # Chimera Graphs
    k_p_list = sys.argv[1].split('_')
    k = k_p_list[1][-1:]
    p = k_p_list[2][:1]

    print(str(k)+"-chimera graph, instance "+str(p))
    print("maxcut(out) : MIQPobj = "+str(miqpobj)+", MIQPcpu ="+str(miqpcpu))
    print("maxcut(out) : MILPobj = "+str(milpobj)+", MILPcpu ="+str(milpcpu))
    print("maxcut(out) : SDPobj = "+str(sdpobj)+", SDPcpu ="+str(sdpcpu))
    print("maxcut(out) : GWobj = "+str(gwobj)+", GWcpu ="+str(gwcpu))
    print('')

    f = open('results_chimera_nonweighted_reduced2k2.out', 'a')
    f.write(str(k)+' '+str(p)+' ')
    f.write(str(miqpobj)+' '+str(milpobj)+' '+str(sdpobj)+'('+str(sdprank)+'/'+str(sdpfull)+') '+str(gwobj)+' ')
    f.write(str(miqpcpu)+' '+str(milpcpu)+' '+str(sdpcpu)+' '+str(gwcpu)+'\n')
    f.close()


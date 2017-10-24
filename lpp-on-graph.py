# Oct 15 2017 Does last passage percolation on lattice.
# Want to compute correlations of busemann functions

import igraph as ig
import numpy as np

#import pdb; pdb.set_trace()

# global variables
# N number of vertices
# g generated igraph
# dbg = 0,1,2,3 debug level.

try:
    dbg
except:
    dbg = 0

randfun = lambda : np.random.exponential(scale=1.0)

def strtuple(s):
    # convert into a tuple
    return tuple( [ int(x) for x in s.split(",") ] )

def tuplestr(*x):
    # flatten allows you to call tuplestr((1,2)) and tuple(1,2)
    global dbg
    x1 = list(flatten(x))
    # turn integers into strings
    x1 = [ str(a) for a in x1 ]
    if dbg >= 3:
        print(x1)
    return ','.join(x1)

def try_addvertex(g,name):
    # will try to add a vertex, if it can't it will fail.
    try:
        g.vs.find(name=name)
    except:
        g.add_vertex(name=name)
    return g

def graphgen(N,directed=True):
    # i want to generate a graph with names. igraph does not check for uniqueness.
    g = ig.Graph(directed=directed)
    g.add_vertex(name=tuplestr(0,0))
    for i in range(0,N):
        for j in range(0,N):
            # i tried to store vertex names as tuples, but it confuses it
            if dbg >= 3 and ((i*j + 1) % 100 == 0):
                print(i,j)
            # try to add two vertices
            g = try_addvertex(g,tuplestr(i+1,j))
            g = try_addvertex(g,tuplestr(i,j+1))
            # once vertices are added if necessary
            g.add_edges([(tuplestr(i,j),tuplestr(i+1,j)),\
                    (tuplestr(i,j),tuplestr(i,j+1))])
    #return ig.Graph.Lattice([N,N],circular=False)
    return g

# finds shortest paths with edge weights. needs 2N(N-1) weights.

# edge weights are easy to use in igraphs library.
# wts = np.random.exponential(size=(2*N*(N-1)))

def vertex_weights(g,wtfun,lpp=True):
    # iterate over vertices, select successor edges for each vertex, and assign edge weight
    # it's a little inefficient
    # must pass a graph and a weight function like np.random.exponential.
    # wtfun must take an optional size argument.
    #Oct 23 2017 Just double the edge weights now, and the sampling algorithm will be 
    # very fast.

    
    wt = list(-wtfun(size=int(g.vcount()/2) ))
    # interleave the two lists to make adjacent edges have equal weights
    wt = [ val for pair in zip(wt) for val in pair ]
    if not lpp:
        # then its fpp; make the weights positive. 
        wt = [ x * (-1) for x in wt ]

    for x in g.vs:
        vwt = wt.pop() # get a weight for a vertex
        for y in g.es.select(_source=x.index):
            y["weight"] = vwt
    return g

def main_loop(number_of_vertices=100):
    global N,g
    N = number_of_vertices
    import pdb; pdb.set_trace()
    g = graphgen(N)
    g = vertex_weights(g)

def find_busemanns(number_of_vertices=100,wtfun=np.random.exponential):
    # takes number of vertices, and weight function. Since this is last passage percolation with positive weights, remember to give a negative weight function. Then one can safely use dijkstra and throw in an extra minus sign to find the last passage time.
    # the extra minus sign is taken into account in the return statement.
    global N,g
    try:
        # hopefully g and N are globally defined. These are used in graph gen.
        N,g
    except:
        g = graphgen(N)
        print("Done generating graph")

    # generate weights
    g = vertex_weights(g,wtfun=wtfun)
    if dbg >= 2:
        print("Done Generating weights")

    # since we're finding last passage times you need to add an extra minus sign in the shortest path.
    times = g.shortest_paths_dijkstra(source=[tuplestr(0,0),tuplestr(1,0),tuplestr(2,0),],target=tuplestr(N-1,N-1),weights='weight')

    # flatten times list using chain
    from itertools import chain
    times = list(chain.from_iterable(times))

    #b2 = g.shortest_paths_dijkstra(tuplestr(2,0),target=tuplestr(N-1,N-1),weights='weight')[0][0] - \
            #g.shortest_paths_dijkstra(tuplestr(1,0),target=tuplestr(N-1,N-1),weights='weight')[0][0]
    # the subtraction takes into account that i've multiplied the weights by -1
    return times[2] - times[1], times[1] - times[0]


def run_find_busemanns(runs=1000, save=True, wtfun=np.random.exponential):
    # runs the find_busemann function several times.
    global bus1,bus2,N,dbg,filename

    bus1 = []
    bus2 = []
    for x in range(0,runs):
        if dbg >= 3 and x % 2 == 0:
            print("run: ",x)
        elif dbg >= 2 and x % 10 == 0:
            print("run: ",x)
        elif dbg >= 1 and x % 50 == 0:
            print("run: ",x)

        p,q = find_busemanns(number_of_vertices=N) 
        bus1.append(p)
        bus2.append(q)

    if save:
        import shelve,datetime
        d = datetime.datetime.today().isoformat()
        filename = 'busemanns-runs-' + str(runs) + '-N-' \
                + str(N) + '-' + d + '.shelf'
        with shelve.open(filename,'c') as shelf:
            shelf['busemanns'] = (bus1,bus2)
            shelf['N'] = N

# busemann functions should have exp(alpha) for horizontal. See romik's lisbook. Recall duality to understand the parameter of the busemann function. E[B] = (\alpha, 1 - \alpha) for \alpha in (0,1). This gives the busemann function with gradient corresponding to -\E[B]. So in the (1,1) direction, one should get the exponential function with parameter 1/2.

def plot_busemann_hist(ret=False):
    global bus1

    h = np.histogram(bus1,density=True)
    # contains the left and right ends of the bins. So h[0] has one less element than h[1]
    x = (h[1][0:len(h[1])-1] + h[1][1:len(h[1])])/2
    y = h[0]
    yexp = np.exp(-x/2)*(1/2)
    l1, = plt.plot(x,y,'-r')
    l2, = plt.plot(x,yexp,'-b')
    plt.legend([l1,l2],['busemann histogram','exp(1/2)'])
    if ret:
        return (x,y,yexp)


def test_indep(bus1,bus2,ind_params=(2.0,2.0),ret=False):
    # indicator funs to test correlations  
    ind12 = lambda x,y: 1.0 if x <= ind_params[0] and y <= ind_params[1] else 0.0
    ind1 = lambda x: 1.0 if x <= ind_params[0] else 0.0
    ind2 = lambda x: 1.0 if x <= ind_params[1] else 0.0

    probs = np.array([ [ind12(a,b),ind1(a),ind2(b)] for a,b in zip(bus1,bus2) ])
    n = len(probs)

    p12 = probs[:,0].sum()/n 
    p1  = probs[:,1].sum()/n 
    p2  = probs[:,2].sum()/n 

    cov = p12  -  p1 * p2
    # i guess variance of indicators is p1 * (1 - p1)
    corr = cov / np.sqrt(p1 * (1-p1) * p2 * (1-p2))

    print('p12,p1,p2,cov,correlation',p12,p1,p2,cov,corr)
    if ret:
        return probs

def import_from_file(filename):
    global bus1,bus2,N

    import shelve
    with shelve.open(filename,'r') as shelf:
        bus1,bus2 = shelf['busemanns']
        try:
            N = shelf['N']
        except:
            print("no N = size of grid defined")
            shelf.close()

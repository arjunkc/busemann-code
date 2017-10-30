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

default_svs = [(0,0),(1,0),(1,0),(2,0)]

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

def vertgen(N):
    i,j = 0,0
    for i in range(0,N):
        for j in range(0,N):
            yield tuplestr(i,j)
    for k in range(0,N):
        yield tuplestr(k,N)
    for k in range(0,N):
        yield tuplestr(N,k)

def edgegen(N):
    i,j = 0,0
    for i in range(0,N):
        for j in range(0,N):
            for k in [(i+1,j),(i,j+1)]:
                yield (tuplestr(i,j),tuplestr(k))

def graphgen(N,directed=True,noigraph_gen=False,asgenerator=True):
    """
    This function creates a directed lattice in d=2 where edges go up or right.
    The ig.Graph.Lattice function does not appear to create directed graphs well.
    Use plot_graph to test with a small N.
    igraph does not check for uniqueness when adding vertices by name.
    Oct 25 2017 The for loop in this function is very slow. An iterator that yields is definitely better
    since the for loop is run by the igraph creation routine.
    Oct 24 2017 This is a fairly inefficient function. Probably easier to add vertices by generating a list of names first.
    noigraph_gen simply retuns edges and vertices
    """

    if asgenerator:
        verts = vertgen(N)
        edges = edgegen(N)
    else:
        verts = []
        edges = []
        for i in range(0,N):
            for j in range(0,N):
                # i tried to store vertex names as tuples, but it confuses it
                if dbg >= 3 and ((i*j + 1) % 100 == 0):
                    print(i,j)
                verts.append(tuplestr(i,j))
                # conditional addition of vertex above i,j
                if i == N-1:
                    verts.append(tuplestr(i+1,j))
                if j == N-1:
                    verts.append(tuplestr(i,j+1))
                edges = edges + [(tuplestr(i,j),tuplestr(i+1,j)),\
                        (tuplestr(i,j),tuplestr(i,j+1))]
    if not noigraph_gen:
        g = ig.Graph(directed=directed)
        g.add_vertices(verts)
        g.add_edges(edges)
        # ig.Graph.Lattice doesn't do directed graphs.
        # return ig.Graph.Lattice([N,N],circular=False)
        return g
    else:
        return verts,edges

def graphgen2(N,directed=True):
    # this is a graph gen version that uses g.add_vertex instead of g.add_vertices.
    g = ig.Graph(directed=directed)
    g.add_vertex(name=tuplestr(0,0))
    for i in range(0,N):
        for j in range(0,N):
            # i tried to store vertex names as tuples, but it confuses it
            if dbg >= 3 and ((i*j + 1) % 100 == 0):
                print(i,j)

            # add vertex to the right of current
            g.add_vertex(tuplestr(i+1,j))
            # since iterating over columns first, do casework to add vertex above it.
            if i == 0:
                g.add_vertex(tuplestr(i,j+1))
            elif j == N-1:
                g.add_vertex(tuplestr(i,j+1))

            # once vertices are added if necessary
            g.add_edges([(tuplestr(i,j),tuplestr(i+1,j)),\
                    (tuplestr(i,j),tuplestr(i,j+1))])

    #return ig.Graph.Lattice([N,N],circular=False)
    return g

def plot_graph(g):
    # testing function that allows you to plot the directed graph
    layout = g.layout_fruchterman_reingold()
    ig.plot(g,layout = layout).show()


def vertex_weights(wtfun,lpp=True):
    # iterate over vertices, select successor edges for each vertex, and assign edge weight
    # could implement this as a generator
    # it's a little inefficient
    # must pass a graph and a weight function like np.random.exponential.
    # wtfun must take an optional size argument.
    #Oct 23 2017 Just double the edge weights now, and the sampling algorithm will be 
    # very fast.

    global N
    # very specific to this graph and the lattice. Can easily be generalized. This particular graph will always have an integer number of vertices.
    # there are N**2 vertices, and each has two outgoing edges.
    wt = list(-wtfun(size=N**2))

    # interleave the two lists to make adjacent edges have equal weights
    wt = [ val for pair in zip(wt,wt) for val in pair ]

    if not lpp:
        # then its fpp; make the weights positive. 
        wt = [ x * (-1) for x in wt ]

    return wt

def find_busemanns(number_of_vertices=100,svs=default_svs,wtfun=np.random.exponential):
    # takes number of vertices, and weight function. Since this is last passage percolation with positive weights, remember to give a negative weight function. Then one can safely use dijkstra and throw in an extra minus sign to find the last passage time.
    # the extra minus sign is taken into account in the return statement.

    # these global variables will be set by this function. should be moved to main_loop
    global N,g,edgewts,mywtfun

    # check graph
    try:
        # hopefully g and N are globally defined. These are used in graph gen.
        g
        # if g exists, check that it came from our lattice generated algorithm with the correct N. The number of edges should be 2N^2 since there are N^2 vertices
        if not 2 * (N ** 2) == g.ecount():
            import time
            print("Started regenerating graph")
            stime = time.time()
            g = graphgen(N)
            print("Ended regenerating graph.",time.time() - stime)
    except:
        import time
        print("Started generating graph.")
        stime = time.time()
        g = graphgen(N)
        print("Ended generating graph.",time.time() - stime)

    # generate weights
    edgewts = vertex_weights(wtfun=wtfun)
    if dbg >= 2:
        print("Done Generating weights")
        print("Generated " + str(len(edgewts)) + " weights")
        print("Using source vertices:", svs)
        print("Using wtfun:", wtfun.__name__)

    # since we're finding last passage times you need to add an extra minus sign in the shortest path.

    times = g.shortest_paths_dijkstra(source=[ tuplestr(x) for x in svs],target=tuplestr(N-1,N-1),weights=edgewts)

    # flatten times list using chain
    from itertools import chain
    times = list(chain.from_iterable(times))
    if dbg >= 2:
        print("times:", times)

    #b2 = g.shortest_paths_dijkstra(tuplestr(2,0),target=tuplestr(N-1,N-1),weights='weight')[0][0] - \
            #g.shortest_paths_dijkstra(tuplestr(1,0),target=tuplestr(N-1,N-1),weights='weight')[0][0]
    # the subtraction takes into account that i've multiplied the weights by -1
    return  times[1] - times[0], times[3] - times[2],


def run_find_busemanns(runs=1000, save=True, number_of_vertices=100, wtfun=np.random.exponential,**kwargs):
    # runs the find_busemann function several times.
    global bus1,bus2,N,dbg,filename,mywtfun,mysvs

    # check N
    try:
        N
    except:
        N = number_of_vertices

    try:
        mywtfun()
        # if I can call this function
    except:
        # set global mywtfun to wtfun
        print("Using weight function:",wtfun.__name__)
        mywtfun = wtfun
    else:
        # then continue
        print("Using weight function:",mywtfun.__name__)

    try:
        mysvs  
    except:
        # if not set
        print("Using default source vertices for busemann functions", default_svs)
        mysvs = default_svs
    else:
        print("Using the following source vertices for busemann functions", mysvs)
    kwargs['svs'] = mysvs

    bus1 = []
    bus2 = []
    import time
    stime = time.time()
    for x in range(0,runs):
        if dbg >= 3 and x % 2 == 0:
            print("run: ",x)
        elif dbg >= 2 and x % 10 == 0:
            print("run: ",x)
        elif dbg >= 1 and x % 50 == 0:
            print("run: ",x)

        p,q = find_busemanns(number_of_vertices=N,wtfun=mywtfun,**kwargs) 
        bus1.append(p)
        bus2.append(q)

    if save:
        save_to_file(runs)

    print("Runtime in seconds: ", time.time() - stime)
    
def print_keys_in_file(f):
    """
    Simple function that prints keys inside a shelf
    Takes 1 argument: a string f containing a filename or a file object.
    """
    with shelve.open(f,'r') as shelf:
        for x in shelf.keys():
            print(x)

def save_to_file(runs):
    global bus1,bus2,mywtfun,N,g,mysvs,default_svs

    import shelve,datetime
    d = datetime.datetime.today().isoformat()

    filename = 'busemanns-runs-' + str(runs) + '-N-' + str(N) \
            + '-' + mywtfun.__name__ + '-' \
            + d + '.shelf'

    with shelve.open(filename,'c') as shelf:
        shelf['busemanns'] = (bus1,bus2)
        shelf['N'] = N
        shelf['wtfun'] = mywtfun.__name__
        shelf['g'] = g
        try:
            mysvs
        except:
            shelf['svs'] = default_svs
        else:
            shelf['svs'] = mysvs


# busemann functions should have exp(alpha) for horizontal. See romik's lisbook. Recall duality to understand the parameter of the busemann function. E[B] = (\alpha, 1 - \alpha) for \alpha in (0,1). This gives the busemann function with gradient corresponding to -\E[B]. So in the (1,1) direction, one should get the exponential function with parameter 1/2.

def plot_busemann_hist(bus1,bins=10,ret=False,bars=False):
    global mywtfun,N


    h = np.histogram(bus1,bins=bins,density=True)
    # contains the left and right ends of the bins. So h[0] has one less element than h[1]
    x = (h[1][0:len(h[1])-1] + h[1][1:len(h[1])])/2
    y = h[0]
    yexp = np.exp(-x/2)*(1/2)

    if ret:
        return (x,y,yexp)

    # otherwise make me a plot
    if not bars:
        l1, = plt.plot(x,y,'-r')
        plt.legend([l1],['busemann histogram, ' + mywtfun.__name__ + ' wts'])
        plt.title('N = ' + str(N))
    else:
        # this density simply normalizes it
        plt.hist(bus1,bins=bins,density=True)
        plt.title('N = ' + str(N) + ', ' + mywtfun.__name__ + ' wts')

    if mywtfun.__name__ == 'exponential':
        l2, = plt.plot(x,yexp,'-b')
        yexp = np.exp(-x/2)*(1/2)
        plt.legend([l2],['exp(1/2)'])

def gammapdf(x,shape,scale):
            return gamma.pdf(x,shape,scale=scale)

def gamma_fit(x,y):
    import scipy.optimize as opt
    from scipy.stats import gamma
    popt,pcov = opt.curve_fit(gammapdf,x,y,p0=(1,1))
    return popt

def test_indep(bus1,bus2,ind_params=(2.0,2.0),ret=False, printout=True):
    # indicator funs to test correlations  
    ind12 = lambda x,y: 1.0 if x >= ind_params[0] and y >= ind_params[1] else 0.0
    ind1 = lambda x: 1.0 if x >= ind_params[0] else 0.0
    ind2 = lambda x: 1.0 if x >= ind_params[1] else 0.0

    probs = np.array([ [ind12(a,b),ind1(a),ind2(b)] for a,b in zip(bus1,bus2) ])
    n = len(probs)

    p12 = probs[:,0].sum()/n 
    p1  = probs[:,1].sum()/n 
    p2  = probs[:,2].sum()/n 

    cov = p12  -  p1 * p2
    # i guess variance of indicators is p1 * (1 - p1)
    if 0 < p1 < 1 and 0 < p1 < 1:
        corr = cov / np.sqrt(p1 * (1-p1) * p2 * (1-p2))
    else:
        corr = nan

    if printout:
        for x,y in zip(['p12','p1','p2','covariance','correlation coeff'],[p12,p1,p2,cov,corr]):
                print(x,': ',y)
    
    if ret:
        return p12,p1,p2,cov,corr

def find_corr(bus1,bus2):
    prd = [ x*y for (x,y) in zip(bus1,bus2) ]
    print("covariance: ", np.mean(prd) - np.mean(bus1)*np.mean(bus2))
    print("correlation coeff: ", np.corrcoef(bus1,bus2))

def plot_correlation(bus1,bus2,plotpoints=20,plottype='covariance',ret=False,printout=False,**kwargs):
    global mywtfun
    minr = max(min(bus1),min(bus2)) 
    maxr = min(max(bus1),max(bus2))
    x = np.linspace(minr,maxr,plotpoints)
    y = []

    # disable debugging
    if dbg >= 2:
        import ipdb
        ipdb.set_trace()

    samples = len(bus1)
    x1 = list(x)
    for a in x:
        p12,p1,p2,cov,corr = test_indep(bus1,bus2,ind_params=(a,a),ret=True,printout=printout)
        if plottype == 'covariance':
            y.append(cov)
        elif plottype == 'correlation':
            if isnan(corr):
                # since x is linspaces, we can assume a is unique
                # this is a dirty rotten hack to remove nans  
                x1.remove(a)
            else:
                y.append(corr)

    l1, = plt.plot(x1,y,**kwargs)
    plt.title(mywtfun.__name__ + ' N=' + str(N) + ', samples=' + str(samples))
    plt.legend([ l1 ],[ plottype ])

    if ret:
        return x1,y

def import_from_file(filename):
    global bus1,bus2,N,mywtfun,g,svs,default_svs,mysvs

    import shelve
    with shelve.open(filename,'r') as shelf:
        bus1,bus2 = shelf['busemanns']
        # number of grid pts
        try:
            N = shelf['N']
        except:
            print("no N = size of grid defined")
            shelf.close()

        # weight function defined. I could have saved the function, or I could have saved the name.
        try:
            shelf['wtfun']
        except:
            print("Error importing wtfun. manually set mywtfun.")
            print("might also need to import numpy if wtfun is from numpy.")
            import traceback
            traceback.print_exc(file=sys.stdout)
        else:
            # there are two types. either callable or just the function name. 
            # this checks for both
            if callable(shelf['wtfun']):
                mywtfun = shelf['wtfun']
            else:
                mywtfun = eval(shelf['wtfun'])
            print("Found wtfun", mywtfun.__name__)

        try:
            mysvs = shelf['svs']
        except:
            mysvs = default_svs
            print("source vertices not set. svs set to default: ", mysvs)
        else:
            print("source vertices found. svs set to: ", mysvs)

        try:
            N = shelf['N']
        except:
            print("no N = size of grid defined")
            shelf.close()
        else:
            print("Found N = ", N)
        # graph
        try:
            shelf['g']
        except:
            print('No graph saved to file')
        else:
            # this is a little slow
            g = shelf['g']
            print('Found graph corresponding to N = ', np.sqrt(g.ecount()/2) )

def absnormal(*args,**kargs):
    return abs(np.random.normal(*args,**kargs))    

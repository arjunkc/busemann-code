# Oct 15 2017 Does last passage percolation on lattice.
# Want to compute correlations of busemann functions

import sys,math,readline

import time
import shelve,datetime

import numpy as np

import igraph as ig

from itertools import chain

import scipy 
import scipy.optimize as opt
from scipy.stats import gamma

import matplotlib.pyplot as plt

# traceback
import traceback

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

def str_to_tuple(s):
    # convert a string like '0,0' into a tuple (0,0)
    return tuple( [ int(x) for x in s.split(",") ] )

def tuple_to_str(*x):
    """
    converts a tuple into a string
    I want it to work for both (0,0) and '0,0' so that its convenient for me
    """
    # flatten allows you to call tuple_to_str((1,2)) and tuple(1,2)
    # flatten doesn't always work in all versions of python. so switched to chain.from_iterable
    global dbg
    # chain.from_iterable, unlike flatten, does not like (1,2). instead it prefers ((1,), (2,))
    # so my list comprehension converts that.
    x1 = list(chain.from_iterable( (i if isinstance(i,tuple) else (i,) for i in x)))
    # turn integers into strings
    x1 = [ str(a) for a in x1 ]
    if dbg >= 3:
        print(x1,','.join(x1))
    return ','.join(x1)

def try_addvertex(g,name):
    # will try to add a vertex, if it can't it will fail.
    try:
        g.vs.find(name=name)
    except:
        g.add_vertex(name=name)
    return g

def vertgen(N):
    """
    generator for graphgen. 
    N = size of N x N 2d square lattice
    will generate vertices from (0,0) to (N-1,N-1)
    """
    i,j = 0,0
    for i in range(0,N):
        for j in range(0,N):
            yield tuple_to_str(i,j)

def edgegen(N):
    """
    generator for graphgen. 
    N = size of N x N 2d square lattice
    will generate 2*(N-1)*N edges  = 
    2 for each vertex (i,j) such that 0<= i < N-1 and 0 <= j < N-1 
    and then 2(N-1) more edges for the two boundary = 2 (N-1)^2 + 2(N-1)
    """
    i,j = 0,0
    for i in range(0,N):
        for j in range(0,N):
            if i != N-1:
                yield (tuple_to_str(i,j),tuple_to_str(i+1,j))
            if j != N-1:
                yield (tuple_to_str(i,j),tuple_to_str(i,j+1))


def graphgen(N,directed=True,noigraph_gen=False,return_object=True):
    """
    This function creates a directed lattice in d=2 where edges go up or right.
    The ig.Graph.Lattice function does not appear to create directed graphs well.
    Use plot_graph to test with a small N.

    igraph does not check for uniqueness when adding vertices by name.

    noigraph_gen does not generate the graph. It works best if asgenerator is set to false so that you get a list of vertices and edges.

    asgenerator=True means verts python generators that "yield" vertices each time it is called. It seemed to me that it might be marginally faster if you passed generators to the igraph.add_vertices(generator) versus igraph.add_vertices( [list of vertices ] ). I haven't really tested this.

    Oct 25 2017 The for loop in this function is very slow. An iterator that yields is definitely better since the for loop is run by the igraph creation routine.

    Oct 24 2017 This is a fairly inefficient function. Probably easier to add vertices by generating a list of names first. noigraph_gen simply retuns edges and vertices
    """

    if dbg >= 3:
        print('running graphgen as generator')
    verts = vertgen(N)
    edges = edgegen(N)

    # make a graph layout for plotting
    if not noigraph_gen:
        if dbg >= 3:
            print('generating new graph')
        try:
            g = ig.Graph(directed=directed)
            g.add_vertices(verts)
            g.add_edges(edges)
        except:
            print('Error in generating igraph')

        # make layout for plotting
        layoutlist = [ (x,y) for x in range(N) for y in range(N) ] 

        if return_object:
            return g,ig.Layout(layoutlist)
        else:
            return g,layoutlist
        #return g,ig.Layout(layoutlist)
    else:
        return ([x for x in verts],[ x for x in edges ])

def plot_graph(g,graphlayout=None,**kwargs):
    """
    testing function that allows you to plot the directed graph
    you can pass it your custom layout. custom layouts can be generated using graphgen
    """
    if graphlayout == None:
        # default layout is a grid.
        graphlayout = g.layout_fruchterman_reingold()
        
    return ig.plot(g,layout = graphlayout,**kwargs)

def vertex_weights(wtfun,N,lpp=True):
    """
    generates vertex weights by making the outgoing edges from each vertex have the same weight
    there are N**2 vertices, and each has two outgoing edges, except for either i=N-1 or j=N-1, in which case we will have only one outgoing edge.
    """
    wt = []
    for i in range(N-1):
        gen = list(-wtfun(size=N-1))
        wt = wt + [ val for pair in zip(gen,gen) for val in pair ]
        # one for the last edge 
        wt = wt + list(-wtfun(size=1))
    # one for the last row of weights
    wt = wt + list(-wtfun(size=N-1))

    if dbg>=3:
        print(len(wt))

    if not lpp:
        # then its fpp; make the weights positive. 
        wt = [ x * (-1) for x in wt ]

    return wt

def find_busemanns(g,N,wtfun,svs,geodesics=False):
    """
    takes number of vertices, and weight function. Since this is last passage percolation with positive weights, remember to give it a negative weight function. Then one can safely use dijkstra and throw in an extra minus sign to find the last passage time.

    the extra minus sign is taken into account in the return statement.

    geodesics = True does not work at the moment since get_all_shortest paths does not take negative weights. We're looking for longest paths here.
    """

    # generate random weights using wtfun
    # returns edge weights made so that all weights incident to an edge are equal
    edgewts = vertex_weights(wtfun,N)

    if dbg >= 2:
        print("Done Generating weights")
        print("Generated " + str(len(edgewts)) + " weights")
        print("Using source vertices:", svs)
        print("Using wtfun:", wtfun.__name__)

    # since we're finding last passage times you need to add an extra minus sign in the shortest path.

    if geodesics == True:
        # this has a slightly different syntax from shortest_paths_dijkstra
        #paths = g.get_shortest_paths(tuple_to_str(N-1,N-1),to=[ tuple_to_str(x) for x in svs],weights=edgewts)
        # get times to all vertices on the lattice
        #targets=g.vs[
        times = g.shortest_paths_dijkstra(source=[ tuple_to_str(x) for x in svs],target=g.vs['name'],weights=edgewts)
    else:
        times = g.shortest_paths_dijkstra(source=[ tuple_to_str(x) for x in svs],target=tuple_to_str(N-1,N-1),weights=edgewts)

    # flatten times list using chain
    times = list(chain.from_iterable(times))
    if dbg >= 2:
        print("times:", times)

    #b2 = g.shortest_paths_dijkstra(tuple_to_str(2,0),target=tuple_to_str(N-1,N-1),weights='weight')[0][0] - \
            #g.shortest_paths_dijkstra(tuple_to_str(1,0),target=tuple_to_str(N-1,N-1),weights='weight')[0][0]
    # the subtraction takes into account that i've multiplied the weights by -1
    if geodesics == True:
        return  times[1] - times[0], times[3] - times[2], paths
    else: 
        return  times[1] - times[0], times[3] - times[2]

def run_find_busemanns(runs=1000, save=True, number_of_vertices=100, wtfun=np.random.exponential,svs=default_svs,geodesics=False):
    """
    Sets the following global variables/parameters for the computation.
    N:  This sets the size of the N by N grid. Can be set using number_of_vertices in run_find_busemanns.
    mysvs: These are different from the default source vertices, also a global. You can set this using the keyword argument svs=[(x,y),(p,q),(a,b),(n,m)]. 
    mywtfun: Set using wtfun= in run_find_busemanns.
    bus1,bus2: Busemann values. Usually set when you import from a file, or you've run a simulation.
    Then it calls find_busemanns runs number of times. Sets global variables first.
    Then it saves to a shelf if save=True
    """
    # set global variables first
    global bus1,bus2,paths,N,dbg,filename,mywtfun,mysvs,g,graphlayout

    # vertices
    N = number_of_vertices

    print("Using weight function:",wtfun.__name__)
    mywtfun = wtfun

    print("Using source vertices for busemann functions", svs)
    mysvs = svs

    # check graph
    try:
        # hopefully g and N are globally defined. These are used in graph gen.
        g
        # if g exists, check that it came from our lattice generated algorithm with the correct N. The number of edges should be 2N^2 since there are N^2 vertices
        if not 2 * (N ** 2) == g.ecount():
            print("Started regenerating graph")
            stime = time.time()
            g,graphlayout = graphgen(N)
            print("Ended regenerating graph.",time.time() - stime)
    except:
        print("Started generating graph.")
        stime = time.time()
        g,graphlayout = graphgen(N)
        print("Ended generating graph.",time.time() - stime)
    else:
        print("Graph found, (hopefully) roughly corresponds to N x N grid with N =  ", int(np.sqrt(g.vcount())))

    # begin running find_busemanns
    bus1 = []
    bus2 = []
    paths = []
    stime = time.time()
    for x in range(0,runs):
        if dbg >= 3 and x % 2 == 0:
            print("run: ",x)
        elif dbg >= 2 and x % 10 == 0:
            print("run: ",x)
        elif dbg >= 1 and x % 50 == 0:
            print("run: ",x)
        elif dbg >= 0 and x % 1000 == 0:
            print("run: ",x)

        if geodesics == True:
            p,q,r = find_busemanns(g,N,wtfun,svs,geodesics=True) 
            paths.append(r)
        else:
            p,q = find_busemanns(g,N,wtfun,svs,geodesics=False) 
        bus1.append(p)
        bus2.append(q)

    if save:
        save_to_file(runs)

    print("Runtime in seconds: ", time.time() - stime)

def return_times(g,N,wtfun,scaled=False,samples=1):
    """
    g: graph
    N: size of grid
    returns occupied vertex indices
    """
    edgewts = vertex_weights(wtfun,N)
    avgtimes = np.zeros(N*N)
    for i in range(samples):
        times = g.shortest_paths_dijkstra(source=tuple_to_str(0,0),target=g.vs['name'],weights=edgewts)
        # flatten times list using chain so that it is a single list
        times = list(chain.from_iterable(times))
        # convert to a list of positive times, 
        # AND also divide by samples to make averaging easier
        times = [ -1 * x / samples for x in times ]
        avgtimes  = avgtimes + times

    # divide by samples to get an average

    return avgtimes

def return_occupied_vertex_coordinates(vertex_list,times,time_threshold,scaled=True,interface=False):
    """
    this is a helper function. returns (scaled) coordinates so it can be plotted easily.
    pass times calculated by return times function
    a lot of the inefficiency could be avoided if graph has vertices ordered by x coordinate. fix graphgen. 
    vertex_list : obtained from g.vs['name']
    times : times for occupying each vertex in vertex list
    time_threshold: a number like 3. returns vertex sites that are smaller than 3.
    interface: returns only the interface between the "wet region" and "dry region"
    """
    #import pdb; pdb.set_trace()

    N = int(np.sqrt(len(vertex_list))) # assuming only N**2 vertices are obtained.
    occupied_vertices = []
    j =  0
    for i in range(N):
        j = 0
        try:
            while j < N and times[i*N + j] <= time_threshold:
                if not interface:
                    z = str_to_tuple(vertex_list[i*N+j])
                    occupied_vertices.append(z)
                j = j + 1
        except:
            print(i,j)
            traceback.print_exc(file=sys.stdout)  

        # append only the last one before j jumped
        if interface and j > 0:
                z = str_to_tuple(vertex_list[i*N+j-1])
                occupied_vertices.append(z)
    # if scaled is true, scale all the vertices
    if scaled:
        occupied_vertices = [(x/time_threshold,y/time_threshold) for x,y in occupied_vertices]

    # if interface=True, occupied vertices will have the interface
    return occupied_vertices

def plot_shape_pyplot(g,wtfun,N,times,compare_with_exponential=True,thresholds=None,interface=False,colors=['red','white'],meansamples=10000,plot_options={'linewidth':2},exp_plot_options={'linestyles':'dashed','linewidth':2}):
    """
    This plots the limit shape 
    times contains first or last passage times to vertices
    colors contains the occupied and unoccupied vertex colors
    returns plots using the igraph library to plot graphs.
    """

    global dbg

    contour_linestyles='dashed'
    contour_linewidth=2

    try:
        samples = wtfun(size=meansamples)
        mean = np.mean(samples)
        std = np.std(samples)
        if dbg>=1:
            print('mean = ',mean)
            print('std = ',std)
        # if no threshold is given, then 
        if thresholds==None:
            # the default threshold works well with scaled since it returns the "limit shape T[nx]
            thresholds = [N * mean/2]

        if dbg>=2:
            print('time threshold:',thresholds)

        vertices = g.vs['name']
        if dbg>=2:
            print('len of vertices',len(vertices))

        plots = []
        for x in thresholds:
            occupied_vertices = return_occupied_vertex_coordinates(vertices,times,x,scaled=True,interface=interface)
            occupied_vertices = np.array(occupied_vertices)
            # color[0] if occupied, color[1] if not
            if interface:
                # if just the interface, might as well use plt.plot
                plots.append(plt.plot(occupied_vertices[:,0],occupied_vertices[:,1],'r-',**plot_options))
            else:
                plots.append(plt.scatter(occupied_vertices[:,0],occupied_vertices[:,1],c=colors[0],s=0.5,edgecolors=None,**plot_options))
    except:
        traceback.print_exc(file=sys.stdout)  

    if compare_with_exponential:
        gridsize = 100
        xaxis = np.linspace(0,1/mean,gridsize)
        yaxis = np.linspace(0,1/mean,gridsize)
        x,y = np.meshgrid(xaxis,yaxis)
        #
        e = np.vectorize(exponential_limit_curve)
        z = e(x,y,mean=mean,std=std)
        plots.append(plt.contour(x,y,z,[1],**exp_plot_options))

    return plots
    #return times,occulib(g,wtfun):

def exponential_limit_curve(x,y,mean=1,std=1):
    """
    pass x and y to it and it returns the universal exponential curve
    see J. Martin 2004 limiting shape for directed percolation models annals of probability.
    see also Martins survey 2005 equation theorem 6.1
    """
    return  mean*(x + y) + 2*std*np.sqrt(x * y)

def list_to_matrix(g,l):
    """
    This is for easier interfacing with igraph
    Converts a list of numbers [ corresponding to each vertex in the graph to 
    to a matrix of numbers m, where m(v[k]) = l[k] 
    """
    vs = g.vs['name']
    N = g.vcount()
    m = np.zeros((N,N))
    for i in range(N):
        v = str_to_tuple(vs[i])
        m[v[0],v[1]] = l[i]

    # return numpy matrix
    return m

def matrix_to_list(g,m):
    """
    Again for easier interfacing with igraph
    convert a matrix to a list 
    Will assume that it is a numpy matrix
    """
    vs = g.vs['name']
    N = int(np.sqrt(g.vcount()))
    l = []
    for i in range(N):
        for j in range(N):
            # uses the fact that the vertices of the graph are arranged in rows in list for
            # (0,0),\ldots,(0,N-1),(1,0),\ldots,(1,N-1)
            l.append(m[i,j])

    return l
        
def plot_geodesics(g,wtfun,layout,N,vcolors=['red','blue','green','gray'],vshapes = ['square','circle','triangle'],svs=[(0,0),(1,0),(2,0)]):
    """
    Will plot geodesics from source vertices specified by source vertices
    Will find times from source vertices to all target vertices.
    Then for each source vertex, you can see which vertices to label
    by following the prescription:  

    Currently only works when default_svs <= 3. Because this is the length of the
    vshapes and vcolors array.
    """
    times_matrices = []
    color_matrices = []
    shape_matrices = []
    plots = []
    edgewts = vertex_weights(wtfun,N)

    default_vcolor = 'white'
    default_shape = 'hidden'
    vsize = int(50.0/np.sqrt(N))
    #plot_options = {'vertex_size':vsize,'edge_arrow_size':0}
    plot_options = {'vertex_size':vsize,'edge_arrow_size':0,'vertex_frame_color':'white'}

    # make a list of times matrices N * N representing a time to each vertex, one for each source vertex
    #import ipdb; ipdb.set_trace()
    for z in range(len(svs)):
        x = svs[z]
        # make an array of the color black
        cs = np.array([ [default_vcolor for i in range(N)] for j in range(N) ])
        # make an array of the vertex shape
        shapes = np.array([ [default_shape for i in range(N)] for j in range(N) ])
        t = g.shortest_paths_dijkstra(source=tuple_to_str(x),target=g.vs['name'],weights=edgewts)
        # flatten and then multiply by 1 to get last passage time
        t = list(chain.from_iterable(t))
        t = [ -1 * x for x in t ]
        tmatrix = list_to_matrix(g,t)
        i = N-1; j = N-1
        # this loop should run in interpreted order N time
        while i != x[0] or j != x[1]:
            cs[i,j] = vcolors[z]
            shapes[i,j] = vshapes[z]
            if i == x[0]:
                j = j - 1
            elif j == x[1]:
                i = i - 1
            elif tmatrix[i-1,j] >= tmatrix[i,j-1]:
                i = i - 1
            else:
                j = j - 1
        #ipdb.set_trace()
        # have to append the last vertex to the list
        cs[x[0],x[1]] = vcolors[z]
        shapes[x[0],x[1]] = vshapes[z]

        # plot graph
        plot = plot_graph(g,layout,vertex_color=matrix_to_list(g,cs),vertex_shape=matrix_to_list(g,shapes),**plot_options)

        times_matrices.append(tmatrix)
        color_matrices.append(cs)
        shape_matrices.append(shapes)
        plots.append(plot)
    # form joint color matrix
    # initialize with last color matrix
    # following is pretty badly optimized code
    joint_color_matrix = cs
    joint_shape_matrix = shapes 
    for i in range(N):
        for j in range(N):
            for z in range(len(svs)):
                if color_matrices[z][i,j] != default_vcolor:
                    # will tend to have the last color that is not the default color
                    joint_color_matrix[i,j] = color_matrices[z][i,j]
                if shape_matrices[z][i,j] != default_shape:
                    # will tend to have the last color that is not the default color
                    joint_shape_matrix[i,j] = shape_matrices[z][i,j]

    joint_color_vector = matrix_to_list(g,joint_color_matrix)
    joint_shape_vector = matrix_to_list(g,joint_shape_matrix)
    plot = plot_graph(g,layout,vertex_color=joint_color_vector,vertex_shape=joint_shape_vector,**plot_options)
    plots.append(plot)

    return plots,joint_color_vector,joint_shape_vector

def plot_shape_igraph(g,layout,wtfun,N,times,thresholds=None,colors=['red','white'],meansamples=10000):
    """
    times contains first or last passage times to vertices
    colors contains the occupied and unoccupied vertex colors
    returns plots using the igraph library to plot graphs.
    looks a little ugly since we do not have as much control as we do in matplotlib.pyplot
    """

    global dbg

    try:
        mean = np.mean(wtfun(size=meansamples))
        # if no threshold is given, then 
        if thresholds==None:
            thresholds = [N/(mean)]

        plots = []
        for x in thresholds:

            if dbg==2:
                print('time threshold:',thresholds)

            vertices = g.vs['name']
            if dbg==2:
                print('len of vertices',len(vertices))

            occupied_vertices_indices = [ i for i in range(len(vertices)) if times[i] <= x ]
            # color[0] if occupied, color[1] if not
            vcolors= [ colors[0]  if i in occupied_vertices_indices else colors[1] for i in range(len(vertices)) ]
            plots.append(plot_graph(g,layout,vertex_color=vcolors))
    except:
        import traceback
        traceback.print_exc(file=sys.stdout)  

    return plots
    #return times,occupied_vertices

def print_keys_in_file(f):
    """
    Simple function that prints keys inside a shelf
    Takes 1 argument: a string f containing a filename or a file object.
    """
    with shelve.open(f,'r') as shelf:
        for x in shelf.keys():
            print(x)

def save_to_file(runs,override_filename=''):
    """
    Saves the following parameters to a file
    g:  graph
    N:  N x N grid size
    busemanns: bus1,bus2 just computed
    wtfun:  The weight function you just used.
    svs:    The source vertices for the Busemann functions.
    """

    global bus1,bus2,mywtfun,N,g,mysvs,default_svs,myfilename

    d = datetime.datetime.today().isoformat()

    # see if global filename set
    if override_filename == '':
        filename = 'busemanns-runs-' + str(runs) + '-N-' + str(N) \
                + '-' + mywtfun.__name__ + '-' \
                + d + '.shelf'
    else:
        filename = override_filename

    with shelve.open(filename,'c') as shelf:
        shelf['busemanns'] = (bus1,bus2)
        shelf['N'] = N
        shelf['wtfun'] = mywtfun.__name__
        try:
            shelf['g'] = g
        except:
            print('Error saving graph. Not saving graph to file.')
        shelf['svs'] = mysvs

def plot_busemann_hist(bus1,bins=10,ret=False,bars=False):
    """
    returns histogram of Busemann function (probability distribution function)
    simply pass a values list to the function and it will plot the pdf using pyplot.hist function
    For exponential vertex weights, busemann functions should have exp(alpha) for horizontal. See romik's lisbook. Recall duality to understand the parameter of the busemann function. E[B] = (\alpha, 1 - \alpha) for \alpha in (0,1). This gives the busemann function with gradient corresponding to -\E[B]. So in the (1,1) direction, one should get the exponential function with parameter 1/2.
    """
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
    """
    These functions were fitted to try and find the pdfs of the Busemann functions
    """
    return gamma.pdf(x,shape,scale=scale)

def gamma_fit(x,y):
    """
    These functions were fitted to try and find the pdfs of the Busemann functions
    """
    popt,pcov = opt.curve_fit(gammapdf,x,y,p0=(1,1))
    return popt

def test_indep(bus1,bus2,ind_params=(2.0,2.0),ret=False, printout=True):
    """
    tests independence, correlation and association of two busemann functions
    """
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
        try:
            corr = cov / np.sqrt(p1 * (1-p1) * p2 * (1-p2))
        except:
            print("nan occurred in test_indep with params: ",ind_params)
    else:
        corr = nan

    if printout:
        for x,y in zip(['p12','p1','p2','covariance','correlation coeff'],[p12,p1,p2,cov,corr]):
                print(x,': ',y)
    
    if ret:
        return p12,p1,p2,cov,corr

def find_corr(bus1,bus2,ret=True):
    prd = [ x*y for (x,y) in zip(bus1,bus2) ]
    cor = np.corrcoef(bus1,bus2)
    cov = np.mean(prd) - np.mean(bus1)*np.mean(bus2)
    print("covariance: ", cov)
    print("correlation coeff: ", cor[0,1])
    if ret:
        return (cov,cor[0,1])

    # commented out code is an attempt to get 3d correlations
    # following does not work as yet.
    #def ind12(x,y,(bus1,bus2)):
        #y = [ 1.0 if b1 >= x and b2 >= y else 0.9 for (b1,b2) in zip(bus1,bus2) ]
        #return sum(y)/len(bus1)

    # following does not work as yet.
    #def ind(x,bus1):
        #y = [ 1.0 if b1 >= x for b1 in bus1 ]
        #return sum(y)/len(bus1)

    #def plot3d_correlation(bus1,bus2,plotpoints=20,plottype='covariance',ret=False,printout=False,**kwargs):
        #global mywtfun
        #minr = max(min(bus1),min(bus2)) 
        #maxr = min(max(bus1),max(bus2))
        #delta = (maxr-minr)/plotpoints
        #xax = np.linspace(minr + delta,maxr - delta,plotpoints)
        #yax = np.linspace(minr + delta,maxr - delta,plotpoints)
        #x,y = np.meshgrid(xax,yax)
    #
        ## lambda and vectorize
        #l12 = lambda x,y : ind12(x,y,(bus1,bus2))
        #e12 = np.vectorize(ind12)
    #
        #l1 = lambda x,y : ind(x,bus1) * ind(y,bus2)
        #e1  = np.vectorize(l1)
    #
        #z12 = e12(x,y) 
        #z1 = e12(x,y)
        #return z12,z1

def plot_correlation(bus1,bus2,plotpoints=20,plottype='covariance',savefig=False,ret=False,printout=False,x_range=None,**kwargs):
    """
    plots correlation for two busemann functions. It plots it for indicators 1{X \leq a} 1{ Y \leq a} and so it only plots something reasonable when X and Y have the same signs for the most part. If X > 0 and Y < 0 this wouldn't work.

    x_range =   tuple (a,b) containing the xrange over which to plot covariances. Useful when you want to compare covariances of two different wtfuns
    """
    global mywtfun
    if x_range == None:
        minr = min(min(bus1),min(bus2)) 
        maxr = max(max(bus1),max(bus2))
        delta = (maxr-minr)/plotpoints
        x = np.linspace(minr + delta,maxr - delta,plotpoints)
    else:
        x = np.linspace(x_range[0],x_range[1],plotpoints)
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
    figtitle = 'busemann ' + plottype + ', ' + mywtfun.__name__ + ' N=' + str(N) + ', samples=' + str(samples) + ', svs=' + str(mysvs)
    plt.title(figtitle)
    #plt.legend([ l1 ],[ plottype ])
    plt.legend() # legend controlled by plt.legend() keyword
    if savefig:
        plt.savefig(figtitle + '.png')

    if ret:
        return l1,x1,y

def import_from_file(filename):
    """
    import saved busemann functions from file
    """
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
            try:
                g = shelf['g']
                print('Found graph corresponding to N = ', np.sqrt(g.ecount()/2) )
            except:
                print('error getting graph edge count')

def absnormal(*args,**kargs):
    """
    absolute value of a normal distribution
    """
    return abs(np.random.normal(*args,**kargs))    

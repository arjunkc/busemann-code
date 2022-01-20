# Oct 15 2017 Does last passage percolation on lattice.
# Want to compute correlations of busemann functions

from functools import wraps
from re import A, S
import sys,math,readline,os

import time
import shelve,datetime

import numpy as np

import igraph as ig

from itertools import chain

import scipy 
import scipy.optimize as opt
from scipy.stats import gamma

import matplotlib.pyplot as plt

# testing
import unittest

# traceback
import traceback

# import cython
# import pyximport; pyximport.install()
# numpy_path = np.get_include()
# os.environ['CFLAGS'] = "-I" + numpy_path
# pyximport.install(setup_args={"include_dirs":numpy_path})
# from compiled import * 
# from compiled import *

exec(open('compiled.pyx').read())

try:
    dbg
except:
    dbg = 0

default_svs = [(0,0),(1,0),(1,0),(2,0)]

def try_addvertex(g,name):
    # will try to add a vertex, if it can't it will fail.
    try:
        g.vs.find(name=name)
    except:
        g.add_vertex(name=name)
    return g

def lpp_num_of_vertices(g,graph_shape='rectangle'):
    """
    returns the size maximal horizontal of a square or triangular grid graph
    That is, if you only have a graph saved, it will produce N. Useful when you load a graph from a file and you want to discover N.
     
    Feb 21 2021 Should be modified to include graph_shape
    """
    # god help you if do not get an integer

    if graph_shape == 'rectangle':
        return math.sqrt(g.vcount())
    elif graph_shape == 'triangle':
        # solve N(N+1)/2 = v
        return (math.sqrt(8*g.vcount() + 1)-1) // 2

def graphgen(N,directed=True,noigraph_gen=False,return_layout_as_object=True,graph_shape='rectangle'):
    """
    This function creates a directed lattice in d=2 where edges go up or right.
    The ig.Graph.Lattice function does not appear to create directed graphs well.
    Use plot_graph to test with a small N.

    Returns: a tuple (g,l) where g is an igraph object, and l is an igraph layout

    directed=True   produces a directed lattice with only up/right paths. Currently the other functions cannot handle the undirected lattice

    noigraph_gen=True does not generate the igraph object. This is mostly used for debugging. 

    return_layout_as_object=True    returns the second return value as an igraph object

    graph_shape='rectangle' or 'triangle'
    If chosen to be a triangle, helps cut down on computation time for limit shape computations. This is because you do not want to limit shape to be truncated.

    igraph does not check for uniqueness when adding vertices by name.

    Oct 25 2017 The for loop in this function is very slow. An iterator that yields is definitely better since the for loop is run by the igraph creation routine.

    Oct 24 2017 This is a fairly inefficient function. Probably easier to add vertices by generating a list of names first. noigraph_gen simply retuns edges and vertices
    """

    if dbg >= 1:
        print('Start generating graph: ' + time.asctime())

    verts = vertgen(N,graph_shape=graph_shape)
    edges = edgegen(N,graph_shape=graph_shape)


    if dbg >= 1:
        print('Done generating vertex and edge lists: ' + time.asctime())

    # make a graph layout for plotting
    if not noigraph_gen:
        if dbg >= 3:
            print('generating new graph')
        try:
            g = ig.Graph(directed=directed)
            #import ipdb; ipdb.set_trace()
            g.add_vertices(verts)
            g.add_edges(edges)
        except:
            print('Error in generating igraph')

        # make layout for plotting
        if graph_shape=='rectangle':
            layoutlist = [ (x,y) for x in range(N) for y in range(N) ] 
        elif graph_shape=='triangle':
            layoutlist = [ (x,y) for x in range(N) for y in range(N - x) ] 

        if dbg >= 1:
            print('Done generating igraph object and layout: ' + time.asctime())

        if return_layout_as_object:
            return g,ig.Layout(layoutlist)
        else:
            return g,layoutlist
    else:
        # if noigraph_gen == True
        return ([x for x in verts],[ x for x in edges ])

def plot_graph(g,graphlayout=None,**kwargs):
    """
    testing function that allows you to plot the directed graph
    you can pass it your custom layout. custom layouts can be generated using graphgen
    """
    if graphlayout == None:
        # default layout is a grid.
        graphlayout = g.layout_fruchterman_reingold()

    width = height = len(g.vs)*20
    kwargs["bbox"] = (width, height)
    kwargs["margin"] = 40
        
    return ig.plot(g,layout = graphlayout,**kwargs)

def find_busemanns(g,N,wtfun,svs,geodesics=False,graph_shape='rectangle'):
    """
    takes number of vertices, and weight function. Since this is last passage percolation with positive weights, remember to give it a negative weight function. Then one can safely use dijkstra and throw in an extra minus sign to find the last passage time.

    the extra minus sign is taken into account in the return statement.

    geodesics = True does not work at the moment since get_all_shortest paths does not take negative weights. We're looking for longest paths here.
    """

    # generate random weights using wtfun
    # returns edge weights made so that all weights incident to an edge are equal
    edgewts = vertex_weights(wtfun,N,graph_shape=graph_shape)

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

def return_times(
        g,
        wtfun=np.random.exponential,
        samples=1,
        graph_shape='rectangle'):
    """
    g: graph
    returns a list of occupied vertex indices
    
    usage: t = return_times(g)

    """
    # import ipdb; ipdb.set_trace() 

    # check if N value correctly set
    try:
        if N != lpp_num_of_vertices(g,graph_shape=graph_shape):
            N = lpp_num_of_vertices(g,graph_shape=graph_shape)
    except:
        N = lpp_num_of_vertices(g,graph_shape=graph_shape)
    
    # Make edge weights by putting the same value of the edge weight on both outgoing edges. This can be slow because of loops in edge generation see vertex_weights()
    if dbg>=1:
        print('Start generating weights: ' + time.asctime())

    # set edge weights. can use wtfun_generator with a wrapper for periodic or vertex weights, or something like np.random.exponential for directly using default edge weights.

    edgewts = -wtfun(size=g.ecount())

    if dbg>=1:
        print('End generating weights: ' + time.asctime())

    avgtimes = np.zeros(g.vcount())  # contains the averaged passage time to each vertex

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

def find_time_threshold(g,N,times):
    """
    The limit shape makes the most sense when you truncate at min(t(0,N-1),t(N-1,0)). However, the times array is one dimensional, so you return the appropriate index

    This now generalizes to the triangular grid
    """

    # find index of leftmost and rightmost vertices
    xmax_vertex = g.vs.find(name='0,'+str(N-1)).index
    ymax_vertex = g.vs.find(name=str(N-1)+',0').index

    return min(times[xmax_vertex],times[ymax_vertex])

def return_occupied_vertex_coordinates(vertex_list,N,times,time_threshold,scaled=True,interface=False,graph_shape='rectangle'):
    """
    this is a helper function. returns (scaled) coordinates so it can be plotted easily.

    Essentially, it returns B_t/t where B_t = { (x,y) : G(x,y) <= t }. Here G is of coarse the passage time. By the LPP analog of the cox durrett theorem, { g(x,y) <= 1 - epsilon } \subset B_t/t \subset { g(x,y) <= 1 + epsilon }.

    pass times calculated by return times function

    a lot of the inefficiency could be avoided if graph has vertices ordered by x coordinate. fix graphgen. 
    vertex_list : obtained from g.vs['name']
    times : times for occupying each vertex in vertex list
    time_threshold: a number like 3. returns vertex sites that are smaller than 3.
    interface: returns only the interface between the "wet region" and "dry region"

    For each x coordinate, you find the largest y such that the occupation time is smaller than the time threshold
    """
    #import pdb; pdb.set_trace()

    occupied_vertices = []
    j =  0
    for i in range(N):
        j = 0
        try:
            # row by row, find the largest vertex smaller than the threshold
            if graph_shape == 'rectangle':
                while j < N and times[i*N + j] <= time_threshold:
                    z = str_to_tuple(vertex_list[i*N+j])
                    if not interface:
                        occupied_vertices.append(z)
                    j = j + 1
            elif graph_shape == 'triangle':
                # num of vertices in i rows
                vert_in_i_rows = i*(2*N-i+1)//2
                while j < N - i and times[vert_in_i_rows + j] <= time_threshold:
                    z = str_to_tuple(vertex_list[vert_in_i_rows+j])
                    if not interface:
                        occupied_vertices.append(z)
                    j = j + 1
        except:
            print(i,j)
            traceback.print_exc(file=sys.stdout)  
        # append only the last vertex if interface
        if interface:
            occupied_vertices.append(z)

    # if scaled is true, scale all the vertices
    if scaled:
        occupied_vertices = [(x/time_threshold,y/time_threshold) for x,y in occupied_vertices]

    # if interface=True, occupied vertices will have the interface
    return occupied_vertices

def plot_shape_pyplot(g,wtfun,N,times,
        compare_with_exponential=True,
        thresholds=None,interface=False,colors=['red','white'],
        meansamples=10000,plot_options={'linewidth':2},
        exp_plot_options={'linestyles':'dashed','linewidth':2},
        graph_shape='rectangle'):
    """
    This plots the limit shape B_t/t where t is chosen to be N * mean/2
    times contains first or last passage times to vertices
    colors contains the occupied and unoccupied vertex colors
    returns plots using the igraph library to plot graphs.
    """

    import ipdb; ipdb.set_trace()

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
            # the default threshold works well with scaled since it returns the limit shape B_t/t. The t used here is the time at which the coordinates (N,0) and (0,N) are hit
            thresholds = [find_time_threshold(g,N,times)]

        if dbg>=2:
            print('time threshold:',thresholds)

        vertices = g.vs['name']
        if dbg>=2:
            print('len of vertices',len(vertices))

        plots = []
        for x in thresholds:
            occupied_vertices = return_occupied_vertex_coordinates(vertices,N,times,x,scaled=True,interface=interface,graph_shape=graph_shape)
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
        # simply plots the limit shape using plt.contour {g_exp(x,y) <= 1}
        gridsize = 100
        xaxis = np.linspace(0,1/mean,gridsize)
        yaxis = np.linspace(0,1/mean,gridsize)
        x,y = np.meshgrid(xaxis,yaxis)
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
    by following the dynamic programming principle:  

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

def plot_shape_igraph(g,layout,wtfun,N,times,
        thresholds=None,colors=['red','white'],
        meansamples=10000):
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

def times_on_diagonal(g,N,times):
    """
    Finds times along antidiagonal line x + y = N - 1
    Uses igraph routines, so works for both triangular and rectangular grid
    """

    verts_diag = [ g.vs.find(name=tuple_to_str(x,N-1-x)).index for x in range(0,N)]
    # find times along diagonal
    times_diag = [ times[x]/N for x in verts_diag ] 
    return times_diag

def plot_exponential_time_constant(
        x,
        wtfun,
        mean=None,
        std=None,
        meansamples=1000000,
        exp_plot_options={'color':'blue','linestyle':'dashed','linewidth':4}):

    """
    plots exponential curve with parameters mean and std on the antidiagonal (x,1-x)

    seems to need at least a million samples to get the mean and standard correct to a few decimal places
    """

    # sample from the weight distribution
    # if both mean and std are not specified
    if mean==None or std==None:
        samples = wtfun(size=meansamples)
        mean = np.mean(samples)
        std = np.std(samples)
        if dbg>=1:
            print('mean = ',mean)
            print('std = ',std)

    #n = plot_points # arbitrary parameter
    #x = np.arange(0,1,1/plot_points)
    e = np.vectorize(exponential_limit_curve)
    y = 1 - x
    z = e(x,y,mean=mean,std=std)
    plt.plot(x,z,**exp_plot_options)
    return x,z

def plot_time_constant(g,wtfun,N,times,plot_options={'linewidth':2,'color':'red'}):

    """
    Simply plots time constant along the diagonal line x + y = N
    """

    times_diag = times_on_diagonal(g,N,times)
        
    x = np.arange(0,1,1/N)
    plt.plot(x,times_diag,**plot_options)
    
    return x,times_diag

def print_keys_in_file(f):
    """
    Simple function that prints keys inside a shelf
    Takes 1 argument: a string f containing a filename or a file object.
    """
    with shelve.open(f,'r') as shelf:
        for x in shelf.keys():
            print(x)

def save_to_file(runs=0,vars_to_save=None,override_filename='',filename_prefix='busemann-runs'):
    """
    Saves the following parameters to a file
    g:  graph

    saves graph g by default, otherwise if vars_to_save is not empty, then it saves

    vars_to_save is a dictionary of the form 
    { 'var_name':var_value, ... }

    automatically appends datetime.shelf to filename

    busemanns, N, wtfun

    busemanns : bus1,bus2 just computed
    N         : N x N grid size
    wtfun     : The weight function you just used.
    svs       : The source vertices for the Busemann functions.
    """

    global bus1,bus2,mywtfun,N,g,mysvs,default_svs,myfilename

    # see if global filename set
    if override_filename == '':
        filename = filename_prefix + str(runs) + '-N-' + str(N) \
                + '-' + mywtfun.__name__  
    else:
        filename = override_filename

    # append datetime and shelf to filename
    d = datetime.datetime.today().isoformat()
    filename = filename + '-' + d + '.shelf'

    with shelve.open(filename,'c') as shelf:
        if vars_to_save == None:
            shelf['busemanns'] = (bus1,bus2)
            shelf['N'] = N
            shelf['wtfun'] = mywtfun.__name__
            shelf['svs'] = mysvs
            try:
                shelf['g']  = g
                # somehow this fails if run this way. I have tried to make it work by saving it in the interpreter, and it does appear to work
            except:
                print('Error saving graph. Not saving graph to file.')
        else:
            for v in vars_to_save:
                shelf[v] = vars_to_save[v]
        shelf.close()

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
    figtitle = plottype + ' ' + mywtfun.__name__ + ' N=' + str(N) + ', samples=' + str(samples) + ' svs=' + str(mysvs)
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

    This has to be done better. This has to read the keys from the shelf and automatically restore it
    """
    global bus1,bus2,N,mywtfun,g,svs,default_svs,mysvs

    import shelve
    with shelve.open(filename,'r') as shelf:
        keys = shelf.keys()
        print("Available keys in shelf are:\n",keys)
        
        if 'busemanns' in keys:
            bus1,bus2 = shelf['busemanns']
            print('Imported busemanns bus1,bus2')
        # number of grid pts
        if 'N' in keys:
            N = shelf['N']
            print('Imported N')

        # weight function defined. I could have saved the function, or I could have saved the name.
        if 'wtfun' in keys:
            # there are two types. either callable or just the function name. 
            # this checks for both
            if callable(shelf['wtfun']):
                mywtfun = shelf['wtfun']
            else:
                mywtfun = eval(shelf['wtfun'])
            print("Found wtfun", mywtfun.__name__)
            print("Set as mywtfun")

        # sometimes there are an array of wtfuns
        if 'wtfuns' in keys:
            wtfuns = shelf['wtfuns']
            print('Imported wtfuns')

        if 'mysvs' in keys:
            mysvs = shelf['svs']

        # graph
        if 'g' in keys:
            g = shelf['g']
            print('Imported graph as g')
            print('Found graph corresponding to N = ', np.sqrt(g.ecount()/2) )

        # times
        if 'times' in keys:
            times = shelf['times']
            print('Imported times as times')

def absnormal(*args,**kargs):
    """
    absolute value of a normal distribution
    """
    return abs(np.random.normal(*args,**kargs))    

def wtfun_generator(g,N,
        periodic_weights=False,
        period=1,
        use_vertex_weights=False,
        set_weight_label_in_graph=True,
        graph_shape='rectangle',
        random_fc = np.random.uniform,
        size=0):
    """ 
    Usage: Define a wrapper, before passing to return_times as follows
    wtfun_wrapper = lambda **x: wtfun_generator(g,N,**x)

    periodic_weights: use periodic weights by repeating a box of weights of size period.
    period:  the size of the period
    vertex_weights: use vertex weights instead of default edge weights
    keyword argument size will not be used

    Jun 15 2021: To do, use the original vertex_weights function, or just incorporate the code from there here.
    """
    if periodic_weights:
        m = period
    else:
        m = N + 1

    if graph_shape == 'rectangle':
        ecount = 2*(N-1)*N
        weights = np.zeros(ecount)

        tempSize = 2*(m-1)*m-2*(m-1)
        tempWeight = random_fc(size=tempSize)
    elif graph_shape == 'triangle':
        weights = np.zeros((N-1)*N)
        sq = 2*(m-1)*m-2*(m-1)
        if use_vertex_weights:
            tempSize = sq // 2
            tempWeight = random_fc(size=tempSize)
        else:
            tempSize = sq
            tempWeight = random_fc(size=tempSize)

    k = 0
    for i in range(m-1):
        for j in range(m-1):
            # i,j -> i+1,j
            if graph_shape == 'rectangle':
                arr = get_idArr(g,i,j,0,m,N)
            elif graph_shape == 'triangle':
                arr = get_idArr(g,i,j,0,m,N,graph_shape='triangle')

            for e in arr:
                weights[e] = tempWeight[k]
            if not use_vertex_weights:
                k = k+1

            # i,j -> i,j+1
            if graph_shape == 'rectangle':
                arr = get_idArr(g,i,j,1,m,N)
            elif graph_shape == 'triangle':
                arr = get_idArr(g,i,j,1,m,N,graph_shape='triangle')
            for e in arr:
                weights[e] = tempWeight[k]
            k = k+1
            
    if set_weight_label_in_graph:
        g.es['label'] = ["{:.3f}".format(weights[i]) for i in range(len(weights))]

    return weights

def  get_idArr(g,i,j,direction,m,N,graph_shape='rectangle'):
    """
    returns an array of edge ids such that all edge weights in this array will share the same weight
    The position of the starting vertices of those edges satisfied x = i+p*(m-1) and y = j+p*(m-1)
    @param direction: distinguish between horizontal and vertical edges, 0 for horinzontal and typically 1 for vertical
    @param m: period
    """
    # initialize array saving eids
    arr = []

    lim = math.ceil((N-1)/(m-1))+1

    if direction == 0: # horizontal
        if graph_shape == 'rectangle':
            xlim = N-1
            ylim = N
        elif graph_shape == 'triangle':
            xlim = N-1-j
            # ylim = N-1-i
        # print(xlim,ylim)
        for p in range(lim):
            for q in range(lim):
                x = i+p*(m-1)
                y = j+q*(m-1)

                if graph_shape == 'triangle':
                    ylim = N-1-x
                # print(xlim,ylim)

                if x < xlim and y < ylim:
                    xName = str(x)+','+str(y)
                    yName = str(x+1)+','+str(y)
                    u = g.vs.find(name=xName).index
                    v = g.vs.find(name=yName).index
                    arr.append(g.get_eid(u,v))
                else:
                    break
    else: # vertical
        if graph_shape == 'rectangle':
            xlim = N
            ylim = N-1
        elif graph_shape == 'triangle':
            # xlim = N-1-i
            ylim = N-1-i

        for p in range(lim):
            for q in range(lim):
                x = i+p*(m-1)
                y = j+q*(m-1)

                if graph_shape == 'triangle':
                    xlim = N-1-y
                
                if x < xlim and y < ylim:
                    xName = str(x)+','+str(y)
                    yName = str(x)+','+str(y+1)
                    u = g.vs.find(name=xName).index
                    v = g.vs.find(name=yName).index
                    arr.append(g.get_eid(u,v))
                else:
                    break
        
    return arr

def gpl(gpp,h): 
    """
    returns gpl given times on the diagonal {(x,y) : x + y = N}
    gpp     :   array containing gpp values diagonal
    N       :   grid size
    h       :   2-tuple containing the "tilt" vector for gpl
    """
    transVerts = [[x/N,(N-1-x)/N] for x in range(0,N)] # scaled vertices on diagonal
    
    # uses the formula gpl = sup( gpp(xi) + h . xi) see equation 4.3 in Georgiou et al
    hp = [np.dot(h,vert) for vert in transVerts]
    pl = np.array(times)+hp
    
    return np.max(pl)

def plot_pl_time_constant(
        gpp,
        hrange=10,
        plotpoints=100,
        **plot_options):
    """
    as the name of the function says, it plots the point to line time constant.
    This is a trivial function, but useful.
    """
        
    x = np.linspace(-hrange,hrange,plotpoints)
    y = [gpl(times,N,[h,-h]) for h in x]

    plt.plot(x,y)

def printA(
        g,
        m,
        arr):
    """
    Designed to print the adjacency matrix of the periodically weighted lattice. We call it A.
    The matrix should have size m^2 x m^2. Each element of the matrix corresponds to a vertex.
    """
    print(' ',end='\t')
    # arr is m^2 x m^2 matrix
    # len(arr) = m^2 of course. Why not just use that?
    for i in range(len(arr)):
        # prints integral part of i/m and the remainder
        # this forms the top row of the array. 
        # it represents coordinates in the lattice graph
        name = str(int(i/m))+','+str(i%m)
        print(name,end='\t')
    print()
    for i in range(len(arr)):
        # the name contains the vertex it is mapped to
        name = str(int(i/m))+','+str(i%m)
        print(name,end='\t')
        for j in range(len(arr)):
            print(format(arr[i][j],'.4f'),end='\t')
        print()

def vertex_to_adjacency_matrix_element(i,j,m):
    """
    The matrix arranges the vertices in the following order: row1,row2,...
    It keeps track of periodicity in the matrix
    i:  horizontal coordinate
    j:  vertical coordinate
    """
    return (i%m)*m+j%m

def adjacency_matrix_element_to_vertex(x,m):
    # The matrix arranges the vertices in the following order: row1,row2,...
    # the x coordinate of the vertex is given by the integer part of x/m
    # the y coordinate of the vertex is given by the remainder of x/m
    return (int(x/m),i%m)

def form_periodic_adj_matrix(g,m):
    """
    Requires a rectangular graph g, with at least m vertices.
    Requires g to have weights assigned to edge labels.

    Forms a periodic adjacency matrix A, by looking at the first m x m box in g,
    and then proceeds to ``periodize'' it.
    """

    num_vertices = m**2

    # form a helper matrix with 0s in it.
    # helper keeps track of which entries of the (sparse) adjacency matrix have an element in it
    helper = np.zeros((num_vertices,num_vertices))
    # initialize the adjacency matrix A to -infty
    A = np.ones((num_vertices,num_vertices))*np.NINF

    for i in range(m):
        for j in range(m):
            # find indices of vertices corresponding to horizontal edges 
            u = g.vs.find(name=str(i)+','+str(j)).index
            v = g.vs.find(name=str(i+1)+','+str(j)).index

            # get edge ids of the edge (u,v)
            eid = g.get_eid(u,v)
            # t will contain the edge weight
            t = float(g.es[eid]['label'])
            
            # The matrix arranges the vertices in the following order: row1,row2,...
            # k1 contains the row of the matrix corresponding to vertex (i,j)
            k1 = vertex_to_adjacency_matrix_element(i,j,m)
            # k2 contains the column of the matrix corresponding to vertex (i+1,j)
            # the edge goes from k1 to k2
            k2 = vertex_to_adjacency_matrix_element(i+1,j,m)
            A[k1][k2] = t
            # helper keeps track of the place where you have non-zero weights in the adjacency matrix, and importantly, whether or not the edge is horizontal or vertical
            helper[k1][k2] = 1

            # Same as above, but for VERTICAL edge
            u = g.vs.find(name=str(i)+','+str(j)).index
            v = g.vs.find(name=str(i)+','+str(j+1)).index

            eid = g.get_eid(u,v)
            t = float(g.es[eid]['label'])
            
            k1 = vertex_to_adjacency_matrix_element(i,j,m)
            k2 = vertex_to_adjacency_matrix_element(i,j+1,m)
            A[k1][k2] = t
            helper[k1][k2] = -1

    return np.array(A),np.array(helper)

def modify_adj_matrix(A,helper,h):
    """
    modify_adj_matrix takes as input a (periodic) adjacency matrix of edge weights.
    
    It modifies the adjacency matrix by adding h to the horizontal weights and -h to the vertical weights.

    @param A:   adjacency matrix
    @param helper: tells you where exactly in the matrix there are non-zero entries
    @param h: new h value
    """

    Anew = A.copy()

    # the helper
    hor = np.where(helper==1)
    ver = np.where(helper==-1)

    for i in range(len(hor[0])):
        u,v = hor[0][i],hor[1][i]
        Anew[u][v] = A[u][v]+h[0]
    for i in range(len(ver[0])):
        u,v = ver[0][i],ver[1][i]
        Anew[u][v] = A[u][v]+h[1]

    return Anew

def add(u,v):
    """
    Does addition in the max-plus algebra.
    This involves taking care of infinities carefully. 
    It uses the convention in Max Plus at Work, Heidergott et al, Ch 1, page 13
    u + v = np.NINF if any one of them is np.NINF

    Does this function need to be here?
    """
    if u == np.NINF or v == np.NINF:
        return np.NINF
    else:
        return u+v

def minus(u,v):
    """
    Does subtraction in the max plus algebra. This is the only tricky convention. According to Heidergott et al, 2.7, page 39, we set np.NINF - np.NINF = 0. This is where the A^* matrix is defined.
    """

    if u == np.NINF and v == np.NINF:
        return 0
    else:
        return u-v

def maxplus_eigenvalue(A):
    """
    Computes max-plus eigenvalue of the square matrix A using Karp's algorithm.

    See page 87, Heidergott et al

    By their convention, A[i,j] = w[j,i]. So we begin by first taking the transpose of the matrix.
    I do not think this makes a huge difference. A and A^T should have the same eigenvalue.

    Only works for irreducible matrices A
    """
    # import ipdb; ipdb.set_trace()
    Atrans = np.transpose(A)
    num_vertices = len(Atrans)
    
    # x is a m x (m+1) vector where m = number of vertices
    x = np.ones((num_vertices,num_vertices+1))*np.NINF
    # # choose arbitrary j∈num_vertices and set x(0) = e_j
    j = np.random.randint(0,num_vertices)
    x[j][0] = 0
    # compute x(k) for k=1,...,num_vertices-1
    for i in range(1,num_vertices+1):
        x[:,i] = maxplus_matrix_dot_vector(Atrans,x[:,i-1])
    # print(x)

    _min = np.zeros(num_vertices)
    for i in range(num_vertices):
        _min[i] = np.min([minus(x[i][num_vertices],x[i][k])/(num_vertices-k) for k in range(num_vertices)])
    return np.max(_min)

class maxplus_eigenvalue_test(unittest.TestCase):
    """
    Tests the maxplus_eigenvalue function on a standard matrix.
    Note that the maxplus algorithm assumes that the A[i,j] = w(j,i). 
    So when we call the testcase, we have to give it the transpose of the matrix from Heidergott.

    """
    def test_on_standard_matrix(self):
        e = np.NINF
        # Heidergott et al, Ex 5.1.1
        A1 = [[e,3,e,1],[2,e,1,e],[1,2,2,e],[e,e,1,e]] 
        A2 = [[1,2,e,7],[e,3,5,e],[e,4,e,3],[e,2,8,e]]
        A3 = [[6,2,e,7],[e,3,5,e],[e,4,e,3],[e,2,8,e]]
        self.assertEqual(maxplus_eigenvalue(np.transpose(A1)),2.5)
        self.assertEqual(maxplus_eigenvalue(np.transpose(A2)),5.5)

def maxplus_matrix_dot_vector(arr,v):
    """
    Matrix times a column vector in the max plus algebra
    Helps in the maxplus_eigenvalue function
    """
    x = np.zeros(len(v))
    for i in range(len(x)):
        x[i] = np.max([add(arr[i][k],v[k]) for k in range(len(v))])

    return x

def plot_gpl_from_range(g,m,hrange,save_figure=False):
    """
    plot gpl in the specified range of h values
    save a figure as well, if the named option is set

    Isn't this function written already?
    """
    A,helper = form_periodic_adj_matrix(g,m)
    x = np.linspace(start=-hrange,stop=hrange,num=1000)

    plt.plot(x,[maxplus_eigenvalue(modify_adj_matrix(A,helper,[h,-h])) for h in x])

    if save_figure:
        plt.savefig('m({})_hrange({}).png'.format(m,hrange))

def perpendicularDistance(x0,x1,u0,u1,v0,v1):
    slope = (v1-u1)/(v0-u0)
    offset = v1-slope*v0

    d = (-slope*x0+x1-offset) / (np.sqrt(slope**2+1))

    return d

def DouglasPeucker(x,arr,epsilon):
    ## Find the point with the maximum distance
    dmax = 0
    index = 0
    end = len(arr)
    for i in range(1,end-1):
        d = perpendicularDistance(x[i],arr[i],x[0],arr[0],x[-1],arr[-1])
        if d < dmax:
            index = i
            dmax = d
    
    ResultList = []
    
    ## If max distance is greater than epsilon, recursively simplify
    if dmax > epsilon:
        ## Recursive call
        recResults1 = DouglasPeucker(x[:index+1],arr[:index+1],epsilon)
        recResults2 = DouglasPeucker(x[index:],arr[index:],epsilon)

        # Build the result list
        ResultList = np.append(recResults1,recResults2,axis=0)
    else:
        ResultList = [[x[0],arr[0]], [x[-1],arr[-1]]]
    
    ## Return the result
    return ResultList

def almost_equal(val1,val2,epsilon):
    """
    checks whether two values are almost equal
    """
    return True if abs(val1-val2) < epsilon else False

def get_num_facets(x,arr,epsilon=1e-4,use_DP=False):
    if use_DP:
        return int(len(DouglasPeucker(x,arr,epsilon))/2)
    else:
        num = 0

        i = 4
        while i < len(arr):
            pre2,pre1,curr,post1,post2 = arr[i-2],arr[i-1],arr[i],arr[i+1],arr[i+2]

            # could be facet 
            if not (almost_equal(pre1,curr,epsilon) and almost_equal(curr,post1,epsilon)):
                if almost_equal(pre1,pre2,epsilon) and almost_equal(post1,post2,epsilon):
                    num += 1
                    i += 2
                elif not almost_equal(post1,post2,epsilon):
                    post3 = arr[i+3]
                    if not almost_equal(post2,post3,epsilon):
                        i += 2
            else:
                i += 1
                        
        return num

def find_bad_index(x,y):
    # import ipdb; ipdb.set_trace()
    bad_index = []

    i = 1
    while i < len(x)-1:
        j = 1
        while i+j < len(x):
            delta1 = (y[i+j]-y[i+j-1])/(x[i+j]-x[i+j-1])
            delta2 = (y[i]-y[i-1])/(x[i]-x[i-1])
            if almost_equal(delta1,delta2,1e-3):
                j += 1
            else:
                break
        if j == 1:
            bad_index.append(i)
        i += j
    
    return bad_index

def to_string(h,gpl):
    delta = dgpl(h,gpl)
        
    for i in range(len(h)):
        if i == 0:
            print('h={}\tgpl={:.5f}'.format(h[i],gpl[i]))
        else:
            print('h={}\tgpl={:.5f}\tdelta={:.5f}'.format(h[i],gpl[i],delta[i-1]))

def fine_tuning(h,gpl,A,helper,m):
    import ipdb; ipdb.set_trace()

    bad_index = find_bad_index(h,gpl)
    while len(bad_index) > 0:
        i = bad_index[0]
 
        curr = h[i]
        new_h = [curr-(h[i]-h[i-1])/2]
        new_eig = [maxplus_eigenvalue(modify_adj_matrix(A,helper,[new_h[0],-new_h[0]]),m)]
        to_string(new_h,new_eig)

        h = h[0:i]+new_h+h[i:]
        gpl = gpl[0:i]+new_eig+gpl[i:]
        to_string(h,gpl)

        bad_index = find_bad_index(h,gpl)
        print(bad_index)

    return h,gpl

def dgpl(h,gpl):
    delta = []
    for i in range(1,len(h)):
        delta.append((gpl[i]-gpl[i-1])/(h[i]-h[i-1]))

    return delta

def plot_comparison_eig_gpl(N,m):
    g, layout = graphgen(N)
    wtfun_wrapper = lambda **x: wtfun_generator(g,N,periodic_weights=True,period=m,**x)
    t = return_times(g,wtfun=wtfun_wrapper) 

    x = [-1+0.2*i for i in range(11)]
    A, helper = form_periodic_adj_matrix(g,m)
    y = [maxplus_eigenvalue(modify_adj_matrix(A,helper,[h,-h])) for h in x] #eig
    y2 = [gpl(times_on_diagonal(g,N,t),N,[h,-h]) for h in x] #gpl

    plt.plot(x,y,label='eig')
    plt.plot(x,y2,label='gpl')
    plt.legend()
    plt.show()

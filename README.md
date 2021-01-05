# About

Runs Dijkstra from the igraph library on the directed two dimensional lattice with weights. 

We're currently testing a hypothesis that Busemann functions ought to be negatively correlated.
Busemann functions are differences of passage times (weighted shortest paths) between 3 vertices (Really it's the infinite version of such differences).

The graph is a $N \times N$ grid. $N$ is a global variable that can be set. You can also set the weight function used for the random weights.

# To run this code

Either 

    import lppsim

Or 

    exec(open('lppsim.py').read())

Originally, the package was designed to compute busemann function distributions, which you invoke as fullows

    bus1,bus2 = run_find_busemanns(runs=1000,number_of_vertices=100)

It computes two pairs of functions bus1 and bus2, which corresponds the following differences of passage times:

$$ 
    B_1 = T(0,Ne) - T(e_1,Ne), \quad B_2 = T(e_1,Ne) - T(2e_1,Ne)
$$

where $e_1 = (1,0)$ and $e = e_1 + e_2$. Here, the passage time $T$ is the maximal sum of weights encountered on up-right paths between the two vertices. The two functions $B_1$ and $B_2$ correspond to Busemann functions. The method will also save the busemann data to a file.

The function

    plot_busemann_histograms()

will plot the density (histogram) of the busemann function. For exponentially distributed weights with rate $1$, this busemann function should be Exp(1/2) distributed. This is because I'm considering passage times from $0$ to $N,N$.

You can test the correlations of the two Busemanns using

    test_indep(bus1,bus2)

This will run it 1000 times on a 100 x 100 grid with exponentially distributed passage times. 

You can also **plot limit shapes** and compare it to the exponential solvable model. To do so,

    N = 1000
    g,layout = graphgen(N)

It is convenient to define N globally before running. If N is defined globally, it will be used, if not `return_times` will set it.

The graph generating algorithm takes a bit of time to generate, especially when you use N >= 1000.

Then you need to specify a **weight function** (wtfun)

    mywtfun = lambda **x: np.random.beta(1.5,2,**x)

By default

    wtfun = np.random.exponential

To get a list of passage times to each vertex, while timing each run:

    import time
    start = time.time()
    t = lppsim.return_times(g,wtfun=mywtfun)
    print(time.time() - start)

To plot the **limit shape** $B_t/t$, where $B_t = \{ T(x,y) \leq t \}$ and $t$ is some dynamically set threshold.

    lppsim.plot_shape_pyplot(g,wtfun,N,t,compare_with_exponential=True,interface=True)

# Useful Functions

1.  `graphgen(N,directed=True,noigraph_gen=False,return_layout_as_object=True)`. generates a graph with a given number of vertices. automatically called by `run_find_busemanns`
1.  `vertex_weights(wtfun,N,lpp=True)`: Generate vertex weights for an NxN nearest neighbor lattice. igraph only allows for edge weights, so this sets the edge $(i,j) \to (i+1,j)$ and $(i,j) \to (i,j+1)$ to be exactly the same. Returns negative weights for last passage percolation so that you can still run a standard first-passage time algorithm that allows for negative weights like Bellman-Ford.
1.  `import_from_file` allows you to pick up things from a shelf with stored date. Two parameters are generally saved: N and the busemann functions.
1.  `plot_graph`: Plots your generated graph. It's a good way to see if a lattice is actually generated for small $N$
1.  `find_busemanns(g,N,wtfun,svs,geodesics=False)`: generates vertex weights using vertgen, finds passage times, then returns adjacent busemann functions for the source vertices list svs. 
1.  `run_find_busemanns(runs=1000, save=True, number_of_vertices=100, wtfun=np.random.exponential,svs=default_svs,geodesics=False)`: A wrapper for `find_busemanns` that sets sensible default weight function, a bunch of samples to average, and so on. Finally it saves to file, and it also tracks and prints runtime. Returns a pair of lists, `[bus1, bus2]` that gives the two busemann increments. The default source vertices are `[(0,0),(1,0),(1,0),(2,0)]`
1.  `return_times(g,wtfun=np.random.exponential,scaled=False,samples=1)`: Another useful function. Returns passage times to the whol grid using a single source `(0,0)`. It takes a samples argument, that can be used to estimate expected passage times instead. 
1.  `plot_shape_pyplot(g,wtfun,N,times,compare_with_exponential=True,thresholds=None,interface=False,colors=['red','white'],meansamples=10000,plot_options={'linewidth':2},exp_plot_options={'linestyles':'dashed','linewidth':2})`: This plots the limit shape $B_t/t$ where $t$ is chosen to be $N * mean/2$ by default. You can plot a growing limit shape by passing a list to thresholds. interface says only plot the interface between the wet and dry region. The `compare_with_exponential` option plots the exponential limit (solvable) shape $\{(x,y) \colon g(x,y) = 1\}$ where $g(x,y) = m (x + y) + 2 \sigma \sqrt{ x y }$. The mean $m$ and $\sigma$ are estimate roughly by sampling from the given `wtfun`. The number of samples is controlled by `meansamples`. 
1.  `plot_geodesics(g,wtfun,layout,N,vcolors=['red','blue','green','gray'],vshapes = ['square','circle','triangle'],svs=[(0,0),(1,0),(2,0)])`: Will plot geodesics from source vertices specified by source vertices Will find times from source vertices to all target vertices.  Then for each source vertex, you can see which vertices to label by following the dynamic programming principle for the last-passage time.
1.  `plot_busemann_hist(bus1,bins=10,ret=False,bars=False)`: Plots a histogram so that you can see the distribution of some quantity. It doesn't have to be a busemann function.

# Notes/Changelog

Jun 21 2019 The problem with `get_all_shortest_paths` is that it requires wtfun to be positive! We have negative weight function since we want to find the LONGEST path. So we're a bit doomed here. So I will have to find the shortest path using the usual backward iteration. If geodesics == True, return all passage times for each of the source vertices. Then move backwards to get the shortest paths.

May 08 2018 Have to implement `mark_groups` to plot shortest paths. I should save a shortest path instance at some point. Or I should write a function called `plot_shortest_path` that takes in a graph, and plots all three shortest paths. You will need to find indices using

    g.vs.find(name=tuplestr)

or you could just use the vertex list returned by the `get_shortest_path` function.

Apr 20 2018 Now `save_to_file` saves files with a .db extension, and when importing, *you have to be careful to drop the extension*. That is, *do not* use filename.shelf.db, and use filename.shelf. There was also a problem of the graph not getting exported since it wasn't declared as a local variable in `run_find_busemanns`.

Apr 20 2018 There was some trouble in the tuplestr function with flatten and `chain.from_iterable`. I had to write code that distinguished between 1,2 as separate arguments and (1,2).

Apr 19 2018 I did some more simulations after talking to Pete, and it seems like negative association only seems to be true for uniform and the source vertices are arranged in a down-right fashion. I saved some image files for the absnormal, and there is a region where its not so clear if things are infact negatively associated.

Feb 24 2018 Could add a global override for the `save_to_file` function. Decided against adding a `set_global_variables` function. There would be a lot of weird hacking for this. I added an argument to `save_to_file` instead. I did some more code cleanup, where I moved global variable setting to the `run_find_busemanns` function.

Feb 24 2018 To do: clean up the precedence of the global variables. It's fine to have N, mysvs, g and mywtfun set as global variables, but if they're passed to the function, one should have the function set these global variables. They should all be handled by `run_find_busemanns` or set global variables for easy debugging. The internally set keyword arguments should have precedence.

Feb 23 2018 mysvs and mywtfun appear to be ugly global variable hacks, but it really helps the shelve function. I don't have to pass all these variables to it. 

Feb 15 2018 Working on plot 3d correlations. The functions ind12 and ind1 do not work, and have been commented out. I'm simply going to test the busemann correlation hypothesis. The busemanns appear to be negatively correlated as far as I can tell.

Nov 01 2017 Plot 3d correlations or contour plot. return to `plot3d_correlations`.

Oct 31 2017 Somehow it's not coming up as positively associated in the simulations. There is a telling correlation coefficient image that indicates that the exponential is independent, more or less. The absnormal and uniform are of course fairly strongly correlated. But at some critical point, the correlation coefficient suddently turns negative. May 09 2018 UPDATE: I dont expect $B(0,e_1)$ and $B(0,e_2)$ to be positively or negatively correlated. Since the arrow at $0$ can "switch", I expect the correlation to be distribution dependent.

Oct 27 2017 Added a bunch of other functions.

    strtuple                string to tuple of the form (i,j) used for coordinates
    tuplestr                tuple (i,j) of to string of the form  "i,j" names used for coordinates
    try_addvertex           to be removed. "tries" to add a vertex to a graph by searching it; is very slow.
    vertgen                 generator containing vertices to be added for 2d directed lattice graph.
    edgegen                 generator containing edges to be added for 2d directed lattice graph.
    graphgen                generates a directed 2d lattice (up-right edges) with N vertices on each axis. Names the  
                            vertices with their coordinates.
    graphgen2               slower graph algorithm that uses g.add_vertex(), but the loop runs in slow pure python.
    plot_graph              testing function that plots the graph generated by graph gen. run with small N graphs.
    vertex_weights          adds vertex weights with given random number generator to graph.
    find_busemanns          find busemann functions given a graph and a random number generator (distribution)
    run_find_busemanns      runs find busemanns runs= number of times and returns a list of busemanns. Calls save to file 
                            and prints diagnostic messages.
    save_to_file            saves busemanns, graph, N, and weight function name used to file
    plot_busemann_hist      plots a histogram for one of the busemanns.
    gammapdf                pdf of the gamma distribution i used for fitting.
    gamma_fit               fit the gamma distribution to the pdf of the busemann function.
    test_indep              test the independence of the two busemann functions; uses 1_{B_1 \geq a} 1_{B_2 \geq a}
    plot_correlation        plots correlation coefficient and covariance of the two busemanns
    import_from_file        imports busemanns, graphs and other things from a 'shelf' file.
    absnormal               absolute value of a normally distributed random variable.

Oct 25 2017 I implemented the vertices and edges as generators, and its an order of a magnitude faster!!

    In [213]: timeit.Timer(stmt='g = graphgen(100)',globals=globals()).timeit(number=1)
    Out[213]: 0.3496230710297823

    In [214]: timeit.Timer(stmt='g = graphgen(100,asgenerator=False)',globals=globals()).timeit(number=1)
    Out[214]: 1.1191907830070704

    In [215]: timeit.Timer(stmt='g = graphgen2(100)',globals=globals()).timeit(number=1)
    Out[215]: 7.132349282968789

Oct 25 2017 Todo: The busemann function for uniform weights looks like a Gamma. Could a gamma with certain parameters fit it?

Oct 25 2017 Todo: check the Busemann function for edge weights. Is it still exponential? Need to write an `edge_weights` function for this. Quite easy to do. Wonder what the Busemann function for uniform is?

Oct 25 2017 Now testing graph generation speeds and seeing what's the slow step. I wonder if it's the igraph step. It certainly seems like it, since graphgen2 took 748 seconds (this calls `g.add_vertex()` everytime) and graphgen that creates the vertex list runs in 176 seconds for N=300. This generally seems to be the case.

    N = 100     graphgen1 1s    graphgen2 6s

The bulk of the time is taken to run through the loop and create the vertex list. So the loop is pretty slow, I guess, maybe because it's pure python. Maybe it's better to get the 

    graph.add_vertices()

algorithm to run through its own loop, since that is probably written in C?

Oct 24 2017 It seems to indicate negative correlations for even the uniform distribution. To check the code to see if it's giving the right answers.

Oct 24 2017 It runs really fast now for about a 1000 runs on a graph of size 100. The step that takes the longest is the graph generation. I'm fine with this for now.

Previously, the slowest step was the vertex weight generation, that took up to a minute. This is because I was doing something stupid.

I should check the paths and see if the paths are sort of independent. If a BK inequality holds, then one should expect that these generally should only align along the geodesic from (1,0) to (N,N).

I must note that the Busemanns rapidly converge to the exponential distribution though! It makes sense, I guess.



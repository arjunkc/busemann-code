# About

Runs Dijkstra from the igraph library on the directed two dimensional lattice with weights. 

We're currently testing a hypothesis that Busemann functions ought to be negatively correlated.
Busemann functions are differences of passage times (weighted shortest paths) between 3 vertices (Really it's the infinite version of such differences).

The graph is a $N \times N$ grid. $N$ is a global variable that can be set. You can also set the weight function used for the random weights.

# To run this code

    exec(open('lpp-on-graph.py').read())
    run_find_busemanns(runs=1000,number_of_vertices=100)

This will run it 1000 times on a 100 x 100 grid with exponentially distributed passage times. 
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

Other useful functions include

1.  graphgen. generates a graph with a given number of vertices. automatically called by `run_find_busemanns`
1.  `vertex_weights`
1.  `import_from_file` allows you to pick up things from a shelf with stored date. Two parameters are generally saved: N and the busemann functions.

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



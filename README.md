# About

Runs Dijkstra from the igraph library on the directed two dimensional lattice with weights. 

We're currently testing a hypothesis that Busemann functions ought to be negatively correlated.
Busemann functions are differences of passage times (weighted shortest paths) between 3 vertices (Really it's the infinite version of such differences).

The graph is a $N \times N$ grid. $N$ is a global variable that can be set. You can also set the weight function using for the random weights.

# To run this code

    exec(open('lpp-on-graph.py').read())
    run_find_busemanns(runs=1000,number_of_vertices=100)

This will run it 1000 times on a 100 x 100 grid with exponentially distributed passage times. 
It computes two pairs of functions bus1 and bus2, which corresponds the following differences of passage times:

$$ B_1 = T(0,Ne) - T(e_1,Ne), \quad B_2 = T(e_1,Ne) - T(2e_1,Ne)$$

where $e_1 = (1,0)$ and $e = e_1 + e_2$. Here, the passage time $T$ is the maximal sum of weights encountered on up-right paths between the two vertices. The two functions $B_1$ and $B_2$ correspond to Busemann functions. The method will also save the busemann data to a file.

The function

    plot_busemann_histograms()

will plot the density (histogram) of the busemann function. For exponentially distributed weights with rate $1$, this busemann function should be Exp(1/2) distributed. This is because I'm considering passage times from $0$ to $N,N$.

You can test the correlations of the two Busemanns using

    test_indep(bus1,bus2)

Other useful functions include

1.  graphgen
1.  `vertex_weights`
1.  `import_from_file` allows you to pick up things from a shelf with stored date. Two parameters are generally saved: N and the busemann functions.

# Notes/Changelog

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



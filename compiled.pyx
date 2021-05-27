#cython: language_level=3

from itertools import chain

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
    # chain.from_iterable, unlike flatten, does not like (1,2). instead it prefers ((1,), (2,))
    # so my list comprehension converts that.
    x1 = list(chain.from_iterable( (i if isinstance(i,tuple) else (i,) for i in x)))
    # turn integers into strings
    x1 = [ str(a) for a in x1 ]

    return ','.join(x1)

def vertex_weights(wtfun,N,graph_shape='rectangle',lpp=True):
    """
    generates vertex weights by making the outgoing edges from each vertex have the same weight.

    Only works on the directed Z^2 lattice, where all edges associated with a source vertex are generated at once

    lpp=True returns -wtfun weights, and lpp=False returns positive wtfun

    Now, you just pass the number of vertices to it, and it will give you a doubled list of weight, which indicates that both weights are equal.

    Originally, I passed N, the size of the rectangular grid. In this case, there are N**2 vertices, and each has two outgoing edges, except for either i=N-1 or j=N-1, in which case we will have only one outgoing edge.
    """
    #import ipdb; ipdb.set_trace()

    if graph_shape == 'rectangle':
        wt = []
        # for loops are very slow here, so you have to make sure you're not running N**2 loop
        # list construction is relatively fast
        for i in range(N-1):
            gen = list(-wtfun(size=N-1))
            wt = wt + [ val for pair in zip(gen,gen) for val in pair ]
            # one for the last edge 
            wt = wt + list(-wtfun(size=1))
        # one for the last row of weights
        wt = wt + list(-wtfun(size=N-1))
    elif graph_shape == 'triangle':
        # generate negative weights
        # num vertices that have outgoing edges = (N-1)*N/2
        num_verts = (N-1)*N // 2
        gen = list(-wtfun(size=num_verts))
        wt = [ val for pair in zip(gen,gen) for val in pair ]

    #if lppsim.dbg>=3:
        #print(len(wt))

    if not lpp:
        # then its fpp; make the weights positive. 
        wt = [ x * (-1) for x in wt ]

    return wt

def vertgen(N,graph_shape='rectangle'):
    """
    generator for graphgen. 
    N = size of N x N 2d square lattice
    will generate vertices from (0,0) to (N-1,N-1)

    graph_shape='rectangle' or 'triangle'
    If chosen to be a triangle, helps cut down on computation time for limit shape computations. This is because you do not want to limit shape to be truncated.
    """
    i,j = 0,0
    if graph_shape == 'rectangle':
        for i in range(0,N):
            for j in range(0,N):
                yield tuple_to_str(i,j)
    elif graph_shape == 'triangle':
        for i in range(0,N):
            for j in range(0,N-i):
                yield tuple_to_str(i,j)

def edgegen(N,graph_shape='rectangle'):
    """
    generator for graphgen. 

    if graph_shape='rectangle'

    if graph_shape='rectangle'
        N = size of N x N 2d square lattice
        will generate 2*(N-1)*N edges  = 
        2 for each vertex (i,j) such that 0<= i < N-1 and 0 <= j < N-1 
        and then 2(N-1) more edges for the two boundary = 2 (N-1)^2 + 2(N-1)

    if graph_shape='triangle'
        will generate all edges for vertices in the triangle with vertices (0,N-1), (N-1,0). All vertices have outgoing edges, except for vertices on the diagonal x + y = N - 1 
    """
    i,j = 0,0
    if graph_shape == 'rectangle':
        for i in range(0,N):
            for j in range(0,N):
                if i != N-1:
                    yield (tuple_to_str(i,j),tuple_to_str(i+1,j))
                if j != N-1:
                    yield (tuple_to_str(i,j),tuple_to_str(i,j+1))
    elif graph_shape == 'triangle':
        for i in range(0,N):
            for j in range(0,N-i):
                if i + j < N-1:
                    # if (i,j) not on the antidiagonal, yield both outgoing edges
                    yield (tuple_to_str(i,j),tuple_to_str(i+1,j))
                    yield (tuple_to_str(i,j),tuple_to_str(i,j+1))

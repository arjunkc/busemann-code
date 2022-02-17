import numpy as np
import unittest

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
        # if this is an NxN graph, this is how many edges it will have.
        # The edges of the box have coordinates (0,0) and (N-1,N-1)
        ecount = 2*(N-1)*N
        weights = np.zeros(ecount)

        tempSize = 2*m**2
        tempWeight = random_fc(size=tempSize)
    elif graph_shape == 'triangle':
        weights = np.zeros((N-1)*N)
        sq = 2*m**2
        if use_vertex_weights:
            tempSize = sq // 2
            tempWeight = random_fc(size=tempSize)
        else:
            tempSize = sq
            tempWeight = random_fc(size=tempSize)

    k = 0
    for i in range(m):
        for j in range(m):
            # i,j -> i+1,j
            if graph_shape == 'rectangle':
                arr = get_idArr(g,i,j,[1,0],m,N)
            # elif graph_shape == 'triangle':
            #     arr = get_idArr(g,i,j,[1,0],m,N,graph_shape='triangle')

            for e in arr:
                weights[e] = tempWeight[k]
            if not use_vertex_weights:
                k = k+1

            # i,j -> i,j+1
            if graph_shape == 'rectangle':
                arr = get_idArr(g,i,j,[0,1],m,N)
            # elif graph_shape == 'triangle':
            #     arr = get_idArr(g,i,j,[0,1],m,N,graph_shape='triangle')
            for e in arr:
                weights[e] = tempWeight[k]
            k = k+1
            
    if set_weight_label_in_graph:
        g.es['label'] = ["{:.3f}".format(weights[i]) for i in range(len(weights))]

    return weights

def is_vertex_in_box(N,i,j,graph_shape='rectangle'):
    if graph_shape == 'rectangle':
        return True if i >= 0 and j >= 0 and i < N and j < N else False
    # elif graph_shape == 'triangle':
    #     return True if i >= 0 and j >=0  and i < N and j < N-1-i else False  
    else:
        print('invalid graph shape')
        return False

def get_idArr(g,i,j,direction,m,N,graph_shape='rectangle'):
    """
    Helped function for periodic weights.

    returns an array of edge ids such that all edge weights in this array will share the same weight

    (i,j) is a vertex in the periodic box.

    The position of the starting vertices of those edges satisfy x = i+p*(m-1) and y = j+p*(m-1)
    @param direction: distinguish between horizontal and vertical edges, 0 for horinzontal and typically 1 for vertical
    @param m: period
    """
    # initialize array saving eids
    arr = []

    # lim has the max number of times the vertex (0,0) repeats in a box with periodicity m, and size N. Ex. N=3,m=2 gives lim=2. The vertex (0,0) repeats twice, but the vertex (0,1) repeats just once.

    lim = math.ceil(N/m)
    for p in range(lim):
        for q in range(lim):
            ux = i+p*m
            uy = j+q*m

            vx,vy = np.add([ux,uy],direction)

            if graph_shape == 'rectangle':
                if is_vertex_in_box(N,ux,uy) and is_vertex_in_box(N,vx,vy):
                    uName = str(ux)+','+str(uy)
                    vName = str(vx)+','+str(vy)
                    u = g.vs.find(name=uName).index
                    v = g.vs.find(name=vName).index
                    arr.append(g.get_eid(u,v))
    
    return arr

def gpl(times,N,h): 
    """
    returns gpl given times on the diagonal {(x,y) : x + y = N}
    gpp     :   array containing gpp values diagonal
    N       :   grid size
    h       :   2-tuple containing the "tilt" vector for gpl
    """
    transVerts = [[x/N,(N-1-x)/N] for x in range(0,N)] # scaled vertices on diagonal
    #transVerts = [[x/N,(N-x)/N] for x in range(0,N)]
    
    # uses the formula gpl = sup( gpp(xi) + h . xi) see equation 4.3 in Georgiou et al
    hp = [np.dot(h,vert) for vert in transVerts]
    pl = np.array(times)+hp
    
    return np.max(pl)

def plot_pl_time_constant(g,N,
        times,
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
        A):
    """
    Designed to print the adjacency matrix of the periodically weighted lattice. We call it A.
    The matrix should have size m^2 x m^2. Each element of the matrix corresponds to a vertex.
    """
    print(' ',end='\t')
    # arr is m^2 x m^2 matrix
    # len(arr) = m^2 of course. Why not just use that?
    for i in range(m**2):
        # prints integral part of i/m and the remainder
        # this forms the top row of the array. 
        # it represents coordinates in the lattice graph
        name = str(int(i/m))+','+str(i%m)
        print(name,end='\t')
    print()
    for i in range(m**2):
        # the name contains the vertex it is mapped to
        name = str(int(i/m))+','+str(i%m)
        print(name,end='\t')
        for j in range(m**2):
            print(format(A[i][j],'.4f'),end='\t')
        print()

def vertex_to_adjacency_matrix_element(t,m):
    """
    The matrix arranges the vertices in the following order: row1,row2,...
    It keeps track of periodicity in the matrix
    t:  t = (i,j)
    """
    return (t[0]%m)*m+t[1]%m

def adjacency_matrix_element_to_vertex(x,m):
    # The matrix arranges the vertices in the following order: row1,row2,...
    # the x coordinate of the vertex is given by the integer part of x/m
    # the y coordinate of the vertex is given by the remainder of x/m
    return (int(x/m),x%m)

def form_periodic_adj_matrix(g,m):
    """
    Requires a rectangular graph g, with at least m vertices.
    Requires g to have weights assigned to edge labels.

    Forms a periodic adjacency matrix A, by looking at the first m x m box in g,
    and then proceeds to ``periodize'' it.
    """

    num_vertices = m**2

    #ipdb.set_trace()

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
            k1 = vertex_to_adjacency_matrix_element((i,j),m)
            # k2 contains the column of the matrix corresponding to vertex (i+1,j)
            # the edge goes from k1 to k2
            k2 = vertex_to_adjacency_matrix_element((i+1,j),m)
            A[k1][k2] = t
            # helper keeps track of the place where you have non-zero weights in the adjacency matrix, and importantly, whether or not the edge is horizontal or vertical
            helper[k1][k2] = 1

            # Same as above, but for VERTICAL edge
            u = g.vs.find(name=str(i)+','+str(j)).index
            v = g.vs.find(name=str(i)+','+str(j+1)).index

            eid = g.get_eid(u,v)
            t = float(g.es[eid]['label'])
            
            k1 = vertex_to_adjacency_matrix_element((i,j),m)
            k2 = vertex_to_adjacency_matrix_element((i,j+1),m)
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

    Does this function need to be here? Jan 20 2022 We can remove this later after discussion with xuchen
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

    def test_on_gpl_eig(self):
        N = 1000
        m = 3
        g, layout = graphgen(N)
        wtfun_wrapper = lambda **x: wtfun_generator(g,N,periodic_weights=True,period=m,**x)
        t = return_times(g,wtfun=wtfun_wrapper) 

        x = [-1+0.2*i for i in range(11)]
        A, helper = form_periodic_adj_matrix(g,3)
        y = [maxplus_eigenvalue(modify_adj_matrix(A,helper,[h,-h])) for h in x] #eig
        y2 = [gpl(times_on_diagonal(g,N,t),N,[h,-h]) for h in x] #gpl

        diff = [np.abs(y[i]-y2[i]) for i in range(len(y))]
        error = 0.01
        if np.max(diff) >= error:
            print('Max difference between gpl (gpp) and eigenvalue is {} and error is {}'.format(np.max(diff),error))

        # plt.figure()
        # plt.title('N={}_m={}_maxdiff={}'.format(N,m,np.max(diff)))
        # plt.plot(x,y,label='eig')
        # plt.plot(x,y2,label='gpl')
        # plt.legend()
        # plt.savefig('periodic gpl versus eigenvalue test N={}_m={}.png'.format(N,m))
        # plt.clf()

def maxplus_matrix_dot_vector(arr,v):
    """
    Matrix times a column vector in the max plus algebra
    Helps in the maxplus_eigenvalue function
    """
    x = np.zeros(len(v))
    for i in range(len(x)):
        x[i] = np.max([add(arr[i][k],v[k]) for k in range(len(v))])

    return x

def plot_gpl_eigenvalue(g,m,hrange,save_figure=False):
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

def plot_periodic_eigenvalue_vs_gpl(N,
        m,
        plotpoints=100,
        g=None
        ):
    """
    plots comparing
    1) periodic lpp gpl computed via forming a large NxN box of mxm weights, and then computing passage times, gpp, and then legendre transforming to find gpl
    2) eigenvalue of the periodic adjacency matrix

    @params
    plotpoints = number of points to plot in comparison plot
    g = pass existing graph to it
    """

    if g == None:
        g, layout = graphgen(N)

    wtfun_wrapper = lambda **x: wtfun_generator(g,N,periodic_weights=True,period=m,**x)
    t = return_times(g,wtfun=wtfun_wrapper) 
    # f = plot_graph(g,layout)
    # f.save('graph.png')

    #x = [-1+0.2*i for i in range(11)]
    x = np.linspace(-1,1,plotpoints)
    A, helper = form_periodic_adj_matrix(g,m)
    y1 = [maxplus_eigenvalue(modify_adj_matrix(A,helper,[h,-h])) for h in x] #eig
    y2 = [gpl(times_on_diagonal(g,N,t),N,[h,-h]) for h in x] #gpl

    diff = [np.abs(y1[i]-y2[i]) for i in range(len(y1))]
    # print(y2)

    plt.figure()
    plt.title('N={}_m={}_maxdiff={}'.format(N,m,np.max(diff)))
    plt.plot(x,y1,label='eig')
    plt.plot(x,y2,label='gpl')
    plt.legend()
    # plt.show()
    plt.savefig(data_dir + '/' + 'periodic_gpl_versus_eigenvalue N={}_m={}.png'.format(N,m))
    #plt.clf()

    try:
        # if layout was generated, return it
        layout
        return g,layout,t,A,x,y1,y2
    except:
        # otherwise g was passed to the function, so no need to return it
        return t,A,x,y1,y2


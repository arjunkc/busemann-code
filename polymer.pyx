#cython: language_level=3

import numpy as np

def new_partition_vector(z,wt_vector):
    z = (np.append([0],z) + np.append(z,0) ) * wt_vector
    #z = (np.append([0],z) + np.append(z,0) ) 
    return z

def return_partition_fn(
        N,
        wtfun=np.random.uniform,
        wts=None,
        wt_type = 'vertex'
        samples=1,
        beta=1,
        ):
    """
    Will return an array of partition functions corresponding to partition functions on the diagonal of a grid.

    So the graph shape is horizontal by default. 

    Use a wrapper for wtfun if you want vertex weights.

    wt_type : 'vertex' or 'edge' to implement
    """

    # number of weights needed = sum of elements in each off diagonal in triangle
    num_wts = np.sum( np.arange(1,N+1) )
        
    try:
        if wts == None:
            # allows you to use the same weight vector for all matrices
            wts = wtfun(size=num_wts*samples)
    except:
        pass

    avgtimes = np.zeros(N) # contains the averaged passage time to the final hypotenuse of the triangular graph 

    #import ipdb; ipdb.set_trace()

    for i in range(samples):
        # initalize partition function
        # removed 1/2
        z = np.array([np.exp(wts[0] * beta )])
        # remove first weight
        wts = wts[1::]
        for j in np.arange(1,N):
            wt_vector = 1/2*np.exp( wts[0:j+1] * beta )
            z = new_partition_vector(z,wt_vector)
            # pop previously used weights
            wts = wts[j+1::]

        # flatten times list using chain so that it is a single list
        avgtimes = avgtimes + z /samples

    # divide by samples to get an average

    return avgtimes



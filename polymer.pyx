#cython: language_level=3

import numpy as np

def new_partition_vector(z,wt_vector,q,scaled=False):
    """
    scaled = True means that it will return a scaled probability vector instead.
    """
    z = (q*np.append([0],z) + (1-q)*np.append(z,0) ) * wt_vector
    #z = (np.append([0],z) + np.append(z,0) ) 
    if scaled:
        # if used as intended, sum(z) will be the ratio of sum(z_n)/sum(z_{n-1})
        return [z/sum(z),sum(z)]
    else:
        return z

def return_partition_fn(
        N,
        wtfun=np.random.uniform,
        wts=None,
        wt_type = 'vertex',
        samples=1,
        full_history=False,
        beta=1,
        q=1/2,
        ):
    """
    Will return an array of partition functions corresponding to partition functions on the diagonal of a grid.

    The graph shape is triangular by default. 

    full_history: return the full triangular array of the partition function
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

    #import ipdb; ipdb.set_trace()
    if full_history:
        avgz = np.zeros(N)
    else:
        avgz = np.zeros(N) # contains the averaged passage time to the final hypotenuse of the triangular graph 
        avglogsumz = 0 #contains log of partition function (sum of z along diagonal)



    for i in range(samples):
        # initalize partition function
        # removed 1/2
        z = np.array([np.exp(wts[0] * beta )])
        logsumz = np.log(sum(z))
        # remove first weight
        wts = wts[1::]
        for j in np.arange(1,N):
            wt_vector = np.exp( wts[0:j+1] * beta )
            # getting scaled partition function 
            z,c = new_partition_vector(z,wt_vector,q,scaled=True)
            logsumz = logsumz + np.log(c)
            # pop previously used weights
            wts = wts[j+1::]

        # flatten times list using chain so that it is a single list
        avgz = avgz + z /samples
        avglogsumz = avglogsumz + logsumz/samples

    # divide by samples to get an average
    return avgz, avglogsumz



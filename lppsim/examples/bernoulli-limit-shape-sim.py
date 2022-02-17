savedir = '/mnt/Core/Research_Work/First_Passage_Percolation/busemann-code/generated-data'

times = []
pvals = [0.3,0.5,0.6,0.7]
for p in pvals:
    import time
    start = time.time()

    #import ipdb; ipdb.set_trace()

    wtfun = lambda **x: np.random.binomial(1,p,**x)
    t = lppsim.return_times(g,wtfun=wtfun)
    print('Time to run' + str(time.time() - start))
    times.append(t)


    fig = plt.figure()
    fig.set_size_inches([15,15/1.33])

    # plot and compare with interface
    plot = lppsim.plot_shape_pyplot(g,wtfun,N,t,
            compare_with_exponential=True,interface=True) 

    # arbitrarily set
    plt.savefig(savedir + 'limit shape bernoulli p=' + str(p) +
            ' versus exponential ' +
            str(N) + 'x' + str(N) + ' grid.svg')


# finally try to save all the things
lppsim.save_to_file(vars_to_save={'times':times,'pvals':pvals},override_filename=savedir + 'limit shape bernoulli times p = ' +
        str(pvals) + ' ' + datetime.datetime.today().strftime('%Y-%m-%dT%H-%M') 
        + '.shelf')


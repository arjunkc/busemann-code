import time,csv,subprocess
start = time.time()

dry_run=False
lppsim.dbg=1
print('Starting at ',time.asctime())
N = 8000
shape='triangle'
g,layout = lppsim.graphgen(N,graph_shape=shape)

times = []
pvals = [0.3,0.5,0.7]
wtfuns = []

for p in pvals:
    #import ipdb; ipdb.set_trace()
    wtfuns.append(('bernoulli p='+str(p),lambda **x: np.random.binomial(1,p,**x)))

# other distributions 
wtfuns = wtfuns + [('uniform',np.random.uniform),
                   ('lognormal',np.random.lognormal),
                   ('exponential',np.random.exponential),
                   ('chisquared k=0.5',lambda **x: np.random.chisquare(0.5,**x)),
                   ('chisquared k=1.5',lambda **x: np.random.chisquare(1.5,**x))
                   ]


# set figure parameters
plt.rc('font',size=16)
# autolayout removes the borders, more appropriate for latex, I think
plt.rc('figure',figsize=[15,15/1.33],autolayout=True)

for wtfun in wtfuns:
    t = lppsim.return_times(g,wtfun=wtfun[1],graph_shape=shape)
    print('Computed times for ' + wtfun[0] + ': ',time.asctime())


    # plot and compare with interface
    if not dry_run:
        plt.figure()
        x,times_diag = lppsim.plot_time_constant(g,wtfun[1],N,t,compare_with_exponential=True) 
        times.append(times_diag)

        # arbitrarily set
        plt.savefig('/mnt/Core/Research_Work/First_Passage_Percolation/busemann-code/time constant ' +
                wtfun[0] + ' versus exponential ' + str(N) + 'x' + str(N) + ' grid.svg')
        plt.savefig('/mnt/Core/Research_Work/First_Passage_Percolation/busemann-code/time constant ' +
                wtfun[0] + ' versus exponential ' + str(N) + 'x' + str(N) + ' grid.pdf')


# finally try to save all the things
if not dry_run:
    lppsim.save_to_file(vars_to_save={'x':x,'times':times,'wtfuns':[x[0] for x in wtfuns]},override_filename='/mnt/Core/Research_Work/First_Passage_Percolation/busemann-code/passage time for various weight functions on main antidiagonal ' + datetime.datetime.today().strftime('%Y-%m-%dT%H-%M') + '.shelf')

    csv_file = '/mnt/Core/Research_Work/First_Passage_Percolation/busemann-code/passage time for various weight functions on main antidiagonal ' + datetime.datetime.today().strftime('%Y-%m-%dT%H-%M') + '.csv'
    f = open(csv_file,'w')
    a=csv.writer(f)
    a.writerow(['x'] + [y[0] for y in wtfuns])
    a.writerows([a for a in zip(*([x] + times))])
    f.close()
    subprocess.call(['lz4','-4',csv_file])
    


# requires lz4 for compression. comment out line if necessary
import lppsim
import time,csv,subprocess
import matplotlib.pyplot as plt
start = time.time()

dry_run=False
lppsim.dbg=1
print('Starting at ',time.asctime())
N = 10
shape='triangle'
save_data = False

try:
    # g or layout may not even be defined
    if lppsim.lpp_num_of_vertices(g,graph_shape=shape) != N:
        print('N and g do not match')
        g,layout = lppsim.graphgen(N,graph_shape=shape)
except: 
    # if N is not the correct number for the graph
    print('N and g do not match')
    g,layout = lppsim.graphgen(N,graph_shape=shape)

savedir = '/mnt/Core/Research_Work/First_Passage_Percolation/busemann-code/generated-data/'

times = []
#pvals = [0.3,0.5,0.7]
pvals = [0.8,0.9]
wtfuns = []

for p in pvals:
    #import ipdb; ipdb.set_trace()
    # note the p=p construct that creates a localizes p, otherwise the array will all contain a lambda with a reference to the last value of p
    wtfuns.append(('bernoulli p='+str(p),lambda p=p,**x: np.random.binomial(1,p,**x)))

# other distributions 
#wtfuns = wtfuns + [('uniform',np.random.uniform),
                   #('lognormal',np.random.lognormal),
                   #('exponential',np.random.exponential),
                   #('chisquared k=0.5',lambda **x: np.random.chisquare(0.5,**x)),
                   #('chisquared k=1.5',lambda **x: np.random.chisquare(1.5,**x))
                   #]


# set figure parameters
plt.rc('font',size=20)
# plt.rc('text',usetex=True) #makes tex like fonts, but takes a while to generate
# autolayout removes the borders, more appropriate for latex, I think
plt.rc('figure',figsize=[15,15/1.33],autolayout=True)
plt.rc('xtick',labelsize=24)
plt.rc('ytick',labelsize=24)

# unclear what these things do
#plt.rc('xtick.minor',top=False)
#plt.rc('xtick.minor',bottom=False)
#plt.rc('ytick.minor',right=False)
#plt.rc('ytick.minor',left=False)

for wtfun in wtfuns:
    t = lppsim.return_times(g,wtfun=wtfun[1],graph_shape=shape)
    print('Computed times for ' + wtfun[0] + ': ',time.asctime())


    # plot and compare with interface
    if not dry_run:
        plt.figure()
        x,times_diag = lppsim.plot_time_constant(g,wtfun[1],N,t) 
        times.append(times_diag)

        # arbitrarily set
        plt.savefig(savedir + 'time constant ' +
                wtfun[0] + ' ' +
                str(N) + 'x' + str(N) + ' grid.svg')
        plt.savefig(savedir + 'time constant ' +
                wtfun[0] + ' ' +
                str(N) + 'x' + str(N) + ' grid.pdf')


# finally try to save all the things
if not dry_run and save_data:
    lppsim.save_to_file(vars_to_save={'x':x,'times':times,'wtfuns':[x[0] for x in wtfuns]},
            override_filename=savedir + 'passage time for various weight functions on main antidiagonal ' +
            datetime.datetime.today().strftime('%Y-%m-%dT%H-%M') + '.shelf')

    csv_file = savedir + 'passage time for various weight functions on main antidiagonal ' + datetime.datetime.today().strftime('%Y-%m-%dT%H-%M') + '.csv'
    f = open(csv_file,'w')
    a=csv.writer(f)
    a.writerow(['x'] + [y[0] for y in wtfuns])
    a.writerows([a for a in zip(*([x] + times))])
    f.close()
    subprocess.call(['lz4','-4',csv_file])
    


import lppsim
import polymer

def mgf(t): return (np.exp(t) - 1)/t
def mgfrw(d): return (np.exp(d) + 1)/2
# there is some weird fudging of N and N-1 because of the way the RW is setup
def fN(d,N): return np.dot([ np.exp(d * x) / ((mgf(1)**N * mgfrw(d))**(N-1)) for x in range(0,N)], polymer.return_partition_fn(N,beta=1)) 
fN(0.1,10)

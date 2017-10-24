# IPython log file

exec(open("lpp-on-graph.py").read())
exec(open("lpp-on-graph.py").read())
x = g.vs[0]
x.successors()
x.index
g.es.indices
x1 = g.es[0]
x1.tuple
g.es.select(_source=0)
xedges = g.es.select(_source=0)
for y in xedges:
    y["weight"] = randfun()
    
g.es.attribute()
g.es.attribute_names()
g.es.attributes
g.es.attributes("weight")
g.es["weight"]
a = lambda x: np.random.exponential
help(np.random.exponential)
a = lambda x: np.random.exponential(scale=1)
a()
a = lambda : np.random.exponential(scale=1)
a()
a()
len(g)
g[0,0]
g[0,2]
len(g.es)
g[0,1
g[0,1]
g[0,0]
g[0,1]
g[0,2]
g[0,3]
g[1,2]
for x in g.vs:
    for y in x.successors():
        print([x.index,y.index])
        
        
x = g.vs[0]
g.es.select(_source=0,_target=1)
a1 = g.es.select(_source=0,_target=1)
g.attributes()
g.Adjacency
g.Adjacency()
g.Adjacency((0,0))
g.Adjacency(0)
g.Adjacency([0])
help(g.Adjacency)
g.vs[0]["name"] = (0,1)
g.attributes()
g.vs.attributes()
g.vs["name"]
y1 = g.es[0]
y1.source
y1.target
exec(open("lpp-on-graph.py").read())
main_loop(number_of_vertices=2)
def tfun(*x):
    type(x)
    
tfun(1,2)
def tfun(*x):
    print(type(x))
    
    
tfun(1,2)
tfun((1,2))
def tfun(*x):
    print(x,type(x))
    
    
tfun(1,2)
tfun((1,2))
def tfun(*x):
    print(x,len(x),type(x))
    
    
tfun((1,2))
def tfun(*x):
    print(x,len(x),type(x))
    return x

    
    
tfun((1,2))
m = tfun((1,2))
m
type(m)
m[0]
m[1]
len(m)
list(m)
flatten(1,2)
help(flatten)
list(flatten(1,2))
help(flatten)
def tfun(*x):
    print(x,len(x),type(x))
    print(flatten(x))
    return x

    
    
m = tfun((1,2))
def tfun(*x):
    print(x,len(x),type(x))
    print(list(flatten(x)))
    return x

    
    
def tfun(*x):
    print(x,len(x),type(x))
    print(list(flatten(x)))
    return x
m = tfun((1,2))
m = tfun(1,2)
exec(open("lpp-on-graph.py").read())
main_loop(number_of_vertices=2)
exec(open("lpp-on-graph.py").read())
main_loop(number_of_vertices=2)
g.vs['name']
','.join([1,2])
exec(open("lpp-on-graph.py").read())
main_loop(number_of_vertices=2)
graphgen()
g.vs['name']
g.vs.indices
exec(open("lpp-on-graph.py").read())
g.vs.indices
main_loop(number_of_vertices=2)
g
g.vs['name']
','.join(['0', '0'])
g.vs[0]['name']
exec(open("lpp-on-graph.py").read())
main_loop(number_of_vertices=2)
exec(open("lpp-on-graph.py").read())
main_loop(number_of_vertices=2)
exec(open("lpp-on-graph.py").read())
main_loop(number_of_vertices=2)
g.vs[0]['name']
g.vs0['name']
g.vs['name']
exec(open("lpp-on-graph.py").read())
main_loop(number_of_vertices=50)
g.es['weight']
len(g.es['weight'])
g.shortest_paths_dijkstra('0,0')
wts = g.shortest_paths_dijkstra('0,0',weights='weight')
times = g.shortest_paths_dijkstra('0,0',weights='weight')
times
len(times)
len(times[0])
times = g.shortest_paths_dijkstra('0,0',target=tuplestr(N,N),weights='weight')
times = g.shortest_paths_dijkstra('0,0',target=tuplestr(N-1,N-1),weights='weight')
times
times = g.shortest_paths_dijkstra('0,0',target=tuplestr(N-2,N-1),weights='weight')
g.shortest_paths_dijkstra('0,0',target=tuplestr(N-2,N-1),weights='weight')
g.vs[tuplestr(48,49)]
g.vs.select(_name=tuplestr(48,49))
g.vs.select(name=tuplestr(48,49))
g.vs.select(name=tuplestr(48,49))['weight']
g.gs.select(_source=tuplestr(48,49))['weight']
g.es.select(_source=tuplestr(48,49))['weight']
g.gs.select(_source=tuplestr(48,49))
g.es.select(_source=tuplestr(48,49))
g.es.select(_source=tuplestr(48,49))[0]
for x in g.es.select(_source=tuplestr(48,49)):
    print(x['weight'])
    
g.es[tuplestr(48,49),tuplestr(49,49))
g.es[tuplestr(48,49),tuplestr(49,49)]
g.es[tuplestr(48,49),tuplestr(49,49)].attributes()
g.es[tuplestr(48,49),tuplestr(49,49)][0]
g.es[tuplestr(48,49),tuplestr(49,49)]['weight']
a = g.es[tuplestr(48,49),tuplestr(49,49)]
dbg = False
a = g.es[tuplestr(48,49),tuplestr(49,49)]
for x in a:
    print(x)
    
a = g.es[tuplestr(48,49),tuplestr(49,49)]
ia
a
a.attributes
a.attributes()
a.get_attribute_values
a.get_attribute_values()
a.get_attribute_values('weight')
g
a = g.es[(tuplestr(48,49),tuplestr(49,49))]
a.
for x in a:
    print(x)
    
a
for x in a:
    print(x.index)
    
a
a.indices
g.es(5,6)
g.es(5)
g.es['weight']
g.es['weight'][-5:-1]
a
g.shortest_paths_dijkstra('0,0',target=tuplestr(N-1,N-1),weights='weight')
g.shortest_paths_dijkstra('0,0',target=tuplestr(N-2,N-1),weights='weight')
g.shortest_paths_dijkstra('0,0',target=tuplestr(N-1,N-2),weights='weight')
exec(open("lpp-on-graph.py").read())
import pdf
import pdb
pdb.run('lpp-on-graph.py')
pdb.run(exec(open('lpp-on-graph.py').read())
)
pdb.run(exec(open("lpp-on-graph.py").read()))
import lpp-on-graph
import lpp-on-graph.py
import 'lpp-on-graph.py'
import '/home/arjun/Core/Research_Work/First_Passage_Percolation/busemann-code/lpp-on-graph.py'
exec(open("lpp-on-graph.py").read())
exec(open("lpp-on-graph.py").read())
exec(open("lpp-on-graph.py").read())
exec(open("lpp-on-graph.py").read())
exec(open("lpp-on-graph.py").read())
exec(open("lpp-on-graph.py").read())
test_correlation()
exec(open("lpp-on-graph.py").read())
test_correlation()
exec(open("lpp-on-graph.py").read())
test_correlation()
test_correlation()
exec(open("lpp-on-graph.py").read())
test_correlation()
test_correlation()
3 % 5
exec(open("lpp-on-graph.py").read())
test_correlation()
exec(open("lpp-on-graph.py").read())
test_correlation()
np.hist(bus1)
matplotlib.pyplot.hist(bus1)
bus1
len(bus1)
exec(open("lpp-on-graph.py").read())
test_correlation()
test_correlation()
exec(open("lpp-on-graph.py").read())
test_correlation()
g
g = ig.Graph.as_directed(g)
g
N = 2
g = graphgen()
g.is_directed()
print(g)
g = Graph.as_directed(g)
g = ig.Graph.as_directed(g)
print(g)
help(ig.Graph.add_edge)
help(ig.Graph.add_edges)
help(ig.Graph.add_edge)
g = ig.Graph()
type(g)
help(g.add_edge)
help(ig)
exec(open("lpp-on-graph.py").read())
test_correlation()
g
print(g)
main_loop()
exec(open("lpp-on-graph.py").read())
main_loop()
g = ig.Graph(directed=True)
N
for i in range(0,N):
    for j in range(0,N):
        # i tried to store vertex names as tuples, but it confuses it
        g.add_vertex(name=tuplestr(i,j))
        if 0 < i:
            g.add_edge(tuplestr(i-1,j),tuplestr(i,j))
        if 0 < j:
            g.add_edge(tuplestr(i,j-1),tuplestr(i,j))
#return ig.Graph.Lattice([N,N],circular=False)
for i in range(0,N):
    for j in range(0,N):
        # i tried to store vertex names as tuples, but it confuses it
        g.add_vertex(name=tuplestr(i,j))
        if 0 < i:
            g.add_edge(tuplestr(i-1,j),tuplestr(i,j))
        if 0 < j:
            g.add_edge(tuplestr(i,j-1),tuplestr(i,j))
#return ig.Graph.Lattice([N,N],circular=False)
N = 2
g = ig.Graph(directed=True)
for i in range(0,N):
    for j in range(0,N):
        # i tried to store vertex names as tuples, but it confuses it
        g.add_vertex(name=tuplestr(i,j))
        if 0 < i:
            g.add_edge(tuplestr(i-1,j),tuplestr(i,j))
        if 0 < j:
            g.add_edge(tuplestr(i,j-1),tuplestr(i,j))
#return ig.Graph.Lattice([N,N],circular=False)
g
g.es
g.es.indices
exec(open("lpp-on-graph.py").read())
exec(open("lpp-on-graph.py").read())
N = 50; g = graphgen()
dbg
N = 100; g = graphgen()
g
dbg = True; N = 100; g = graphgen()
dbg = 1; N = 100; g = graphgen()
exec(open("lpp-on-graph.py").read())
dbg = 1; N = 100; g = graphgen()
50 % 100
exec(open("lpp-on-graph.py").read())
dbg = 1; N = 100; g = graphgen()
exec(open("lpp-on-graph.py").read())
dbg = 1; N = 100; g = graphgen()
dbg = 1; N = 50; g = graphgen()
dbg = 1; N = 60; g = graphgen()
dbg = 1; N = 80; g = graphgen()
dbg = 1; N = 90; g = graphgen()
dbg = 0
exec(open("lpp-on-graph.py").read())
test_correlation(runs=1)
bus1
g.es.indices
exec(open("lpp-on-graph.py").read())
g
del g
test_correlation()
dbg = 1
del g
test_correlation()
del g
exec(open("lpp-on-graph.py").read())
N=2; test_correlation()
N=2; test_correlation()
exec(open("lpp-on-graph.py").read())
N=2; test_correlation()
print(g)
exec(open("lpp-on-graph.py").read())
N=3; test_correlation()
del g
N=3; test_correlation()
g.es.indices
del g; N=3; test_correlation()
del g; N=3; test_correlation(runs=1)
print(g)
bus1
bus2
del g; N=100; dbg = 2; test_correlation(runs=10)
del g; N=10; dbg = 2; test_correlation(runs=10)
del g; N=50; dbg = 2; test_correlation(runs=10)
del g; N=50; dbg = 1; test_correlation(runs=10)
del g; N=50; dbg = 1; test_correlation(runs=1000)
get_ipython().magic('ls -al')
len(bus1)
plt.hist(bus1)
h = np.histogram(bus1,bins=10,density=True)
h
plt.plot(h)
h
plt.plot(h[0],h[1])
h
plt.plot(h[1],h[0])
h[0]
h[1]
len(h[1])
len(h[0])
h = np.histogram(bus1,density=True)
h
len(h[0])
len(h[0],h[1][0:len(h[1])-1)
len(h[0],h[1][0:len(h[1])-1)])
h[0],h[1][0:len(h[1])-1)]
h[0],h[1][0:len(h[1])-1]
plt.plot(h[0],h[1][0:len(h[1])-1])
x = h[0]
y = h[1][0:len(h[1])-1]
y = h[0]
x = h[1][0:len(h[1])-1]
yexp = np.exp(x)
yexp
yexp = np.exp(-x/2) (1/2)
yexp = np.exp(-x/2)*(1/2)
yexp
plt.plot(x,y,x,yexp)
plt.plot(x,y,'-r',x,yexp)
plt.plot(x,y,'-r',x,yexp,'-b')
get_ipython().magic('logstart')
x = h[1][0:len(h[1])-1]
x = (h[1][0:len(h[1])-1] + h[1][1:len(h[1])])/2
yexp = np.exp(-x/2)*(1/2)
plt.plot(x,y,'-r',x,yexp,'-b')
plt.plot(x,y,'-r',x,yexp,'-b')
plt.legend()
plt.plot(x,y,'-r',label='histogram',x,yexp,'-b',label='exp(1/2)')
plt.plot(x,y,'-r',label=['histogram','exp(1/2)'],x,yexp,'-b')
plt.plot(x,y,'-r',x,yexp,'-b',label=['histogram','exp(1/2)'])
plt.legend()
plt.plot(x,y,'-r',x,yexp,'-b',label='histogram,exp(1/2)')
plt.legend()
plt.plot(x,y,'-r',x,yexp,'-b',label='histogram',label='exp(1/2)')
def plot_busemann_hist(ret=False):
    global bus1

    h = np.histogram(bus1,density=True)
    # contains the left and right ends of the bins. So h[0] has one less element than h[1]
    x = (h[1][0:len(h[1])-1] + h[1][1:len(h[1])])/2
    y = h[1]
    yexp = np.exp(-x/2)*(1/2)
    l1 = plt.plot(x,y,'-r')
    l2 = plt.plot(x,yexp,'-b')
    plt.legend([l1,l2],['busemann histogram','exp(1/2)'])
    plt.legend()
    if ret:
        return (x,y,yexp)
plot_busemann_hist()
def plot_busemann_hist(ret=False):
    global bus1

    h = np.histogram(bus1,density=True)
    # contains the left and right ends of the bins. So h[0] has one less element than h[1]
    x = (h[1][0:len(h[1])-1] + h[1][1:len(h[1])])/2
    y = h[1]
    yexp = np.exp(-x/2)*(1/2)
    l1 = plt.plot(x,y,'-r')
    l2 = plt.plot(x,yexp,'-b')
    plt.legend([l1,l2],['busemann histogram','exp(1/2)'])
    plt.legend()
    if ret:
        return (x,y,yexp)
def plot_busemann_hist(ret=False):
    global bus1

    h = np.histogram(bus1,density=True)
    # contains the left and right ends of the bins. So h[0] has one less element than h[1]
    x = (h[1][0:len(h[1])-1] + h[1][1:len(h[1])])/2
    y = h[0]
    yexp = np.exp(-x/2)*(1/2)
    l1 = plt.plot(x,y,'-r')
    l2 = plt.plot(x,yexp,'-b')
    plt.legend([l1,l2],['busemann histogram','exp(1/2)'])
    plt.legend()
    if ret:
        return (x,y,yexp)
plot_busemann_hist()
bus3 = bus1
bus4 = bus1
import shelve
filename = 'busemanns-runs-' + str(runs) + '-N-' + str(N) + '.shelf'
with shelve.open(filename,'c') as shelf:
    shelf['busemanns'] = (bus1,bus2)
import shelve
filename = 'busemanns-runs-' + str(1000) + '-N-' + str(N) + '.shelf'
with shelve.open(filename,'c') as shelf:
    shelf['busemanns'] = (bus1,bus2)
with shelve.open('shelf-example', 'r') as shelf:
 for key in shelf.keys():
     print(repr(key))  
with shelve.open(filename, 'r') as shelf:
 for key in shelf.keys():
     print(repr(key))  
shelf =  shelve.open(filename, 'r')
shelf.keys
bus3,bus4 = shelf['busemanns']
len(bus3)
bus3[0:10]
ind12 = lambda x,y: 1.0 if x <= params[0] and y <= params[1] else 0.0
params
params = (2.0,2.0)
exec(open("lpp-on-graph.py").read())
test_indep(bus1,bus2)
test_indep(bus1,bus2)
exec(open("lpp-on-graph.py").read())
test_indep(bus1,bus2)
exec(open("lpp-on-graph.py").read())
test_indep(bus1,bus2)
exec(open("lpp-on-graph.py").read())
test_indep(bus1,bus2)
exec(open("lpp-on-graph.py").read())
probs = test_indep(bus1,bus2,ret=True)
probs[0,:]
exec(open("lpp-on-graph.py").read())
probs = test_indep(bus1,bus2,ret=True)
type(g)
times = g.shortest_paths_dijkstra(source=[tuplestr(0,0),tuplestr(1,0),tuplestr(2,0),],target=tuplestr(N-1,N-1),weights='weight')
times
flatten(times)
import datetime
datetime.date.today().isoformat()
datetime.datetime.today().isoformat()
from iterable import chain
from itertools import chain
chain.from_iterable(times)
list(chain.from_iterable(times))
exec(open("lpp-on-graph.py").read())
test_correlation(runs=2)
plot_busemann_hist()
(x,y,yexp) = plot_busemann_hist()
x
y
bus1
bus2
np.histogram(bus1,density=True)
(x,y,yexp) = plot_busemann_hist(ret=True)
x
y
yexp
test_indep(bus1,bus2)
import_from_file('busemanns-runs-1000-N-50.shelf')
exec(open("lpp-on-graph.py").read())
exec(open("lpp-on-graph.py").read())
import_from_file('busemanns-runs-1000-N-50.shelf')
with shelve.open(filename, 'w') as shelf:
     shelf['N'] = N
   
filename
with shelve.open('busemanns-runs-1000-N-50.shelf', 'w') as shelf:
     shelf['N'] = N
   
shelf 
f = close('busemanns-runs-1000-N-50.shelf')
f = open('busemanns-runs-1000-N-50.shelf')
with shelve.open(f, 'w') as shelf:
    shelf['N'] = N
      
f.close()
with shelve.open('busemanns-runs-1000-N-50.shelf', 'w') as shelf:
     shelf['N'] = N
   
with shelve.open('busemanns-runs-1000-N-50.shelf', 'w') as shelf:
     for key in shelf:
         print(key)
     
   
with shelve.open('busemanns-runs-1000-N-50.shelf', 'w') as shelf:
    for key in shelf:
        print(key)
        
      
with shelve.open('busemanns-runs-1000-N-50.shelf', 'r') as shelf:
    for key in shelf:
        print(key)
        
      
with shelve.open('busemanns-runs-1000-N-50.shelf') as shelf:
    for key in shelf:
        print(key)
        
      
shelf 
      
shelf.close()

      
shelf =  shelve.open(filename, 'w')
shelf['N'] = N
shelf.close()
shelf =  shelve.open('busemanns-runs-1000-N-50.shelf', 'w')
shelf =  shelve.open('busemanns-runs-1000-N-50.shelf', 'w')
shelf =  shelve.open('busemanns-runs-1000-N-50.shelf') 
f
f.close()
shelf =  shelve.open('busemanns-runs-1000-N-50.shelf') 
f.close()
f.close()
get_ipython().magic('clear (f)')
dir()
shelf
shelf.close()
shelf
shelf.dict
shelf.clear
shelf.items
shelf.items()
for key in shelf:
    print(key)
    
quit()

# extract a column of x and y coordinates
x_verts = [ lppsim.str_to_tuple(x)[0] for x in g.vs['name'] ]
y_verts = [ lppsim.str_to_tuple(x)[1] for x in g.vs['name'] ]

# make an array of rows to write to a csv file
to_write = [ [x,y,t1] for x,y,t1 in zip(x_verts,y_verts,t) ]

f = open('occupation times ' + mywtfun.__name__ + datetime.date.today().isoformat() + ' 5000 x 5000 grid.csv','w')

# write to csv
import csv
a = csv.writer(f)
a.writerows(to_write)
f.close()


from numpy import savetxt,column_stack,arange

'''
This module generates 2D voronoi diagrams from a list of points.
'''

def voronoi2D(a,P,L):
    # save the data in the right format     
    savetxt(
    'cpp/mpmdata',
    column_stack((arange(len(a)), a)),
    fmt='%d %g %g %g'
    )
 
    # trigger the shell command
    import subprocess
    p = subprocess.Popen(['cpp/import',str(P.G.x_m),str(P.G.x_M),str(P.G.y_m),str(P.G.y_M),'-0.5','0.5'])
    p.wait()
     
    # open the results file and parse results
    f = open('cpp/mpmvoro', 'r')
    for line in f:
        l = line.rstrip().split(' ')
        ID = int(l[0])
        sides = int(l[1])
        L.S[ID].voronoi_volume = l[2]
        L.S[ID].areas = l[3:3+sides]
        L.S[ID].neighbours = l[3+sides:]
        #print l
        #print L.S[ID].voronoi_volume, L.S[ID].neighbours, L.S[ID].areas

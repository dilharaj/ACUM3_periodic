import numpy as np
from array import array

print("\n\n\n\n###### PRE-PROCESSING SCRIPT FOR ACUM-HAMSTR ######\n\n\nThis is the generalized version of the pre-processing script to convert HAMSTR surfall*.dat solutions to ACUM-readable surface inputs\n\nA surface_*.dat file will be created for each body\n\n")

sol_first = int(input("Enter first solution number: "))
sol_last = int(input("Enter last solution number: "))
nbody_first = int(input("Enter first body number: "))
nbody_last = int(input("Enter last body number: "))


num_sols = int(sol_last - sol_first + 1)
nbody =  int(nbody_last - nbody_first + 1)

print("Number of solutions = %d\nNumber of surfaces = %d\n\n"%(num_sols,nbody))

for i in range(nbody):
    
    body_i = str(nbody_first + i)
    N = np.zeros((nbody,num_sols))   # Number of nodes
    E = np.zeros((nbody,num_sols)) # Number of elements

    print("\n############################\n\n        Surface %s\n"%(body_i))
    
    for j in range(num_sols):
        sol_i = str(sol_first + j).zfill(3)
        file = open("surfall"+sol_i+"_"+body_i+".dat")
            
        lst = []

        for line in file:
            lst += [line.split()]

        file.close()

        N[i][j] = int(lst[1][3])
        E[i][j] = int(lst[1][5])
        
        n = int(N[i][j])
        e = int(E[i][j])

        if (j==0):
            print('Number of Nodes =  %d \nNumber of elements = %d\n'%(N[i][0],E[i][0]))
            print('Solution number ...\n')
        if (N[i][j] != N[i][0] or E[i][j] != E[i][0]):
            print('Solution number %s of Body number %s is inconsistent ...! Terminating the process' % (sol_i, body_i))

        print(j)

        x =  np.zeros(n)
        y =  np.zeros(n)
        z =  np.zeros(n)
        
        if (j==0):
            u =  np.zeros((num_sols,e))
            v =  np.zeros((num_sols,e))
            w =  np.zeros((num_sols,e))
            cp = np.zeros((num_sols,e))
            rho = np.zeros((num_sols,e))

        conn = np.zeros((e,4))
        for k in range(n):
            ix = k+3
            iy = k+3+n+1
            iz = k+3+2*(n+1)

#            if (j==16 and k==47633):
#                print(iz)
            x[k] = float(lst[ix][0])
            y[k] = float(lst[iy][0])
            z[k] = float(lst[iz][0])


        for l in range(e):
            irho = l+3+3*(n+1)
            iu = irho+e
            iv = iu+e
            iw = iv+e
            icp = iw+e
            iconn = icp+2*e+1
            
            rho[j][l] = float(lst[irho][0])
            u[j][l] = float(lst[iu][0])
            v[j][l] = float(lst[iv][0])
            w[j][l] = float(lst[iw][0])
            cp[j][l] = float(lst[icp][0])
            
            for q in range(4):
                conn[l][q] = int(lst[iconn][q])

        #print(x[0],y[0],z[0],u[0],v[0],w[0],cp[0],conn[0])
        
        if (j==0):
            X =  np.zeros((num_sols,e))
            Y =  np.zeros((num_sols,e))
            Z =  np.zeros((num_sols,e))
            A =  np.zeros((num_sols,e))
            Nx = np.zeros((num_sols,e))
            Ny = np.zeros((num_sols,e))
            Nz = np.zeros((num_sols,e))
            
        ## Element Area computation
        if (j==0):
            totalarea = 0.0

        for h in range(e):
            x1 = x[int(conn[h][0]-1)]
            y1 = y[int(conn[h][0]-1)]
            z1 = z[int(conn[h][0]-1)]
            x2 = x[int(conn[h][1]-1)]
            y2 = y[int(conn[h][1]-1)]
            z2 = z[int(conn[h][1]-1)]
            x3 = x[int(conn[h][2]-1)]
            y3 = y[int(conn[h][2]-1)]
            z3 = z[int(conn[h][2]-1)]
            x4 = x[int(conn[h][3]-1)]
            y4 = y[int(conn[h][3]-1)]
            z4 = z[int(conn[h][3]-1)]

            s1 = np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2) 
            s2 = np.sqrt((x3-x2)**2+(y3-y2)**2+(z3-z2)**2)
            s3 = np.sqrt((x4-x3)**2+(y4-y3)**2+(z4-z3)**2)
            s4 = np.sqrt((x1-x4)**2+(y1-y4)**2+(z1-z4)**2)
            A[j][h] = (s1*s2+s3*s4)/2
    
            if (j==0):
                totalarea = totalarea + A[j][h]
            
            v1x = x3-x1
            v1y = y3-y1
            v1z = z3-z1
            v2x = x4-x2
            v2y = y4-y2
            v2z = z4-z2

            n1 = v1y*v2z - v1z*v2y
            n2 = v1z*v2x - v1x*v2z
            n3 = v1x*v2y - v1y*v2x
            nmag = np.sqrt(n1**2 + n2**2 + n3**2)
            s = 0.5*nmag

#            if(j==2 and h==10):
#                print(s,A[j][h])

            Nx[j][h] = n1/nmag
            Ny[j][h] = n2/nmag
            Nz[j][h] = n3/nmag

        # Element coordinate computation
        for r in range(e):
            X[j][r] = 0.25*(x[int(conn[r][0]-1)] + x[int(conn[r][1]-1)] + x[int(conn[r][2]-1)] + x[int(conn[r][3]-1)])
            Y[j][r] = 0.25*(y[int(conn[r][0]-1)] + y[int(conn[r][1]-1)] + y[int(conn[r][2]-1)] + y[int(conn[r][3]-1)]) 
            Z[j][r] = 0.25*(z[int(conn[r][0]-1)] + z[int(conn[r][1]-1)] + z[int(conn[r][2]-1)] + z[int(conn[r][3]-1)]) 
        

        #print(X[0],Y[0],Z[0])

           
    print('\nReading writing data for body number %s...\n'%(body_i))


    file = open("surface_"+body_i+".dat",'w')
    
    file.write("%d %d\n" % (e, num_sols))
    for jj in range(num_sols):
        for ii in range(e):
            file.write("%e %e %e %e %e %e %e %e\n" % (X[jj][ii], Y[jj][ii], Z[jj][ii], cp[jj][ii], A[jj][ii], Nx[jj][ii], Ny[jj][ii], Nz[jj][ii]))
                
    file.close

print('Done..!')

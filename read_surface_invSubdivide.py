import numpy as np
from array import array
from scipy.interpolate import UnivariateSpline
import math
import sys
import random
from random import seed
from random import randint
def mean(list):
    return sum(list)/len(list)

def rms(list):
    mean1 = sum(list)/len(list)
    n=0
    for i in range(len(list)):
        n=n+(list[i]-mean1)**2

    return np.sqrt(n/len(list))


print("\n\n\n\n###### PRE-PROCESSING SCRIPT FOR ACUM-HAMSTR ######\n\n\nThis is the generalized version of the pre-processing script to convert HAMSTR surfall*.dat solutions to ACUM-readable surface inputs\n\nA surface_*.dat file will be created for each body\n\n")


print("@ @ @ # # Elements are inverse-subdivided using emap_*.dat files # # @ @ @\n\n")

sol_first = int(input("Enter first solution number: "))
sol_last = int(input("Enter last solution number: "))
nbody_first = int(input("Enter first body number: "))
nbody_last = int(input("Enter last body number: "))


num_sols = int(sol_last - sol_first + 1)
nbody =  int(nbody_last - nbody_first + 1)
E_check = np.zeros((nbody+1,2))


print("Number of solutions = %d\nNumber of surfaces = %d\n\n"%(num_sols,nbody))

for i in range(nbody):
    
    body_i = str(nbody_first + i)
    N = np.zeros((nbody,num_sols))   # Number of nodes
    E = np.zeros((nbody,num_sols)) # Number of elements

    print("\n############################\n\n        Surface %s\n"%(body_i))
    
    for j in range(1):
        sol_i = str(sol_first + j).zfill(3)
        file = open("surfall_time"+sol_i+"_domain"+body_i+".dat")
            
        lst = []

        for line in file:
            lst += [line.split()]

        file.close()

        N[i][j] = int(lst[1][3])
        E[i][j] = int(lst[1][5])
        
        n = int(N[i][j])
        e = int(E[i][j])


    emap_file = "emap_"+str(nbody_first+i-1)+".dat"
    if path.exists(emap_file):
        file = open("emap_"+str(nbody_first+i-1)+".dat",'r')
        lst_map = []
        for line in file:
            lst_map += [line.split()]

        file.close
    else:
        print("No emap file found")
        quit()

    maplen = 0
    for j in range(e):
        if(int(lst_map[j][0])>maplen):
            maplen = int(lst_map[j][0])

    maplen += 1
    print("Number of Clustes : %d\n"%(maplen)) 
    if(len(lst_map) != e):
        print("Mismatch of element numbers...!\n\nExiting...")

        exit()

    emap = np.zeros((maplen,5))
    cmap = np.zeros((maplen,16))
    conn_map = np.zeros((maplen,4))   
   
    #extract connectivity data for all cells
    for j in range(e):
        jj = int(lst_map[j][0])
        emap[jj][0] += 1
        emap[jj][int(emap[jj][0])] = int(j)


    for j in range(num_sols):
        sol_i = str(sol_first + j).zfill(3)
        file = open("surfall_time"+sol_i+"_domain"+body_i+".dat")
            
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
            cp_mean = np.zeros((e))
            cp_rms = np.zeros((e))
            rho = np.zeros((num_sols,e))
            
            uc =  np.zeros((num_sols,maplen))
            vc =  np.zeros((num_sols,maplen))
            wc =  np.zeros((num_sols,maplen))
            cpc = np.zeros((num_sols,maplen))
            rhoc = np.zeros((num_sols,maplen))
        


        conn = np.zeros((e,4))
        for k in range(n):
            ix = k+3
            iy = k+3+n+1
            iz = k+3+2*(n+1)

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
                

        if (j==0):
            X =  np.zeros((num_sols,e))
            Y =  np.zeros((num_sols,e))
            Z =  np.zeros((num_sols,e))
            A =  np.zeros((num_sols,e))
            Nx = np.zeros((num_sols,e))
            Ny = np.zeros((num_sols,e))
            Nz = np.zeros((num_sols,e))

            Xc =  np.zeros((num_sols,maplen))
            Yc =  np.zeros((num_sols,maplen))
            Zc =  np.zeros((num_sols,maplen))
            Ac =  np.zeros((num_sols,maplen))
            Nxc = np.zeros((num_sols,maplen))
            Nyc = np.zeros((num_sols,maplen))
            Nzc = np.zeros((num_sols,maplen))
            

            
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

            Nx[j][h] = n1/nmag
            Ny[j][h] = n2/nmag
            Nz[j][h] = n3/nmag

        # Element coordinate computation
        for r in range(e):
            X[j][r] = 0.25*(x[int(conn[r][0]-1)] + x[int(conn[r][1]-1)] + x[int(conn[r][2]-1)] + x[int(conn[r][3]-1)])
            Y[j][r] = 0.25*(y[int(conn[r][0]-1)] + y[int(conn[r][1]-1)] + y[int(conn[r][2]-1)] + y[int(conn[r][3]-1)]) 
            Z[j][r] = 0.25*(z[int(conn[r][0]-1)] + z[int(conn[r][1]-1)] + z[int(conn[r][2]-1)] + z[int(conn[r][3]-1)]) 
        

        count = 0
        for icc in range(maplen):
            for iel in range(int(emap[icc][0])):

                Xc[j][icc] += X[j][int(emap[icc][iel+1])]
                Yc[j][icc] += Y[j][int(emap[icc][iel+1])]
                Zc[j][icc] += Z[j][int(emap[icc][iel+1])]
                Nxc[j][icc] += Nx[j][int(emap[icc][iel+1])]
                Nyc[j][icc] += Ny[j][int(emap[icc][iel+1])]
                Nzc[j][icc] += Nz[j][int(emap[icc][iel+1])]
                Ac[j][icc] += A[j][int(emap[icc][iel+1])]
                uc[j][icc] += u[j][int(emap[icc][iel+1])]
                vc[j][icc] += v[j][int(emap[icc][iel+1])]
                wc[j][icc] += w[j][int(emap[icc][iel+1])]
                rhoc[j][icc] += rho[j][int(emap[icc][iel+1])]
                cpc[j][icc] += cp[j][int(emap[icc][iel+1])]
                

              

            Xc[j][icc] = Xc[j][icc]/emap[icc][0]
            Yc[j][icc] = Yc[j][icc]/emap[icc][0]
            Zc[j][icc] = Zc[j][icc]/emap[icc][0]
            Nxc[j][icc] = Nxc[j][icc]/emap[icc][0]
            Nyc[j][icc] = Nyc[j][icc]/emap[icc][0]
            Nzc[j][icc] = Nzc[j][icc]/emap[icc][0]
            uc[j][icc] = uc[j][icc]/emap[icc][0]
            vc[j][icc] = vc[j][icc]/emap[icc][0]
            wc[j][icc] = wc[j][icc]/emap[icc][0]
            rhoc[j][icc] = rhoc[j][icc]/emap[icc][0]
            cpc[j][icc] = cpc[j][icc]/emap[icc][0]
            
            Nxc[j][icc] = Nxc[j][icc]/np.sqrt(Nxc[j][icc]**2+Nyc[j][icc]**2+Nzc[j][icc]**2)
            Nyc[j][icc] = Nyc[j][icc]/np.sqrt(Nxc[j][icc]**2+Nyc[j][icc]**2+Nzc[j][icc]**2)
            Nzc[j][icc] = Nzc[j][icc]/np.sqrt(Nxc[j][icc]**2+Nyc[j][icc]**2+Nzc[j][icc]**2)
           
          
    print('\nReading data for body number %s...\n'%(body_i))

     # compute Cp_mean and Cp_rms

#    for r in range(e):
#        cp_mean[r] = mean(cp[:,r])
#        
#        cp_rms[r] = rms(cp[:,r])
#
#    for r in range(maplen):
#        r2 = r+1
#        randnum = random.randrange(100)
#        for rr in range(int(coalist[r2][0])):
#            cp_mean[int(coalist[r2][rr+1])-1] = randnum

    file = open("surface_"+body_i+".dat",'w')
   
    
    file.write("%d %d\n" % (maplen, num_sols))
    for jj in range(num_sols):
        for ii in range(maplen):
             
            file.write("%e %e %e %e %e %e %e %e\n" % (Xc[jj][ii], Yc[jj][ii], Zc[jj][ii], cpc[jj][ii], Ac[0][ii], Nxc[0][ii], Nyc[jj][ii], Nzc[jj][ii]))

    file.close



print('Done..!')

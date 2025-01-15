import numpy as np
import os
import array as array
import math
import shutil


f1 = [0,0,0,1,2,3,3,3,3,0,0,0,0,0,0,0,0,15,16,17,17,17,17,0,0,0,0,0]

n = len(f1)

j = 0
for i in range(n):
	print(i)	
	if f1[i] == 0:
		j += 1
		file = open('../../surfaces/surface_'+str(i+1)+'.dat')
		lst = []
		for line in file:
			lst += [line.split()]
		file.close()
		nel0 = int(lst[0][0])
		iel = nel0
		ntime = int(lst[0][1])
		data = np.zeros((ntime,nel0*5,8))
		for ii in range(ntime):
		    for jj in range(nel0):
			for kk in range(8):
			    data[ii][jj][kk] = float(lst[1+ii*nel0+jj][kk])

		for k in range(n):
		    if f1[k] == i+1:
		        				
			file = open('../../surfaces/surface_'+str(k+1)+'.dat')
			lst = []
			for line in file:
				lst += [line.split()]
			file.close()
			nel1 = int(lst[0][0])

			if(ntime != int(lst[0][1])):
			   print("Inconsistent # of time steps!")
			   sys.exit()
			
			for ii in range(ntime):
			    for jj in range(nel1):
				for kk in range(8):
				    data[ii][iel+jj][kk] = float(lst[1+ii*nel1+jj][kk])

			iel += nel1

		file = open('surface_'+str(j)+'.dat','w')
		file.write("%d %d\n"%(iel,ntime))
		for ii in range(ntime):
		    for jj in range(iel):
			for kk in range(8):
			    file.write("%e "%(data[ii][jj][kk]))
			file.write("\n")

		file.close()





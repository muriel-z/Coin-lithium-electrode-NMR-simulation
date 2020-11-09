#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 07:22:30 2020

@author: muri
"""


import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D 
import calculateFieldShift as cFS 
from SimulationVolume import *

#Voy a construir una celda cilindrica, con las dimensiones h=0.71mm
# y d=12mm. Además quiero FOVz=10*h y FOVxy=2*d


N = [256,256,256]
# Sólo por las dimensiones pedidas y el nro elegido de N, tomo este voxelsize
voxelSize = [0.03, 0.1, 0.1]# mm
volumen = SimulationVolume(voxelSize=voxelSize, N=N)
matriz_3D = np.zeros(N)


#Ordeno las variables necesarias

vsz, vsy, vsx = voxelSize
Nmz, Nmy, Nmx = N
d = 12 #mm
h = 0.71 #mm
R = d/2
ind_R = int(R/vsx)
nsx = int(d/vsx)
nsy = int(d/vsy)
nsz = int(h/vsz)

# Defino los indices que conforman la celda 

indices = []

ind_y = 0
ind_x = 0
iz = 0
 
for ind_y in range(Nmy):
     for ind_x in range(Nmx):
         if (ind_x-Nmx/2)**2 + (ind_y-Nmy/2)**2 < ind_R**2:
             for iz in range(int(Nmz/2)-int(nsz/2), int(Nmz/2)+int(nsz/2)+1):
                 indices.append((iz, ind_y, ind_x))

#%%
# Una vez que tenemos los indices, ponemos un 1 en esos lugares en la matriz
    
indices = np.array(indices).T
indices_flat = np.ravel_multi_index(indices, N)   
np.put(matriz_3D, indices_flat, 1)
#%%
# ahora calculamos la perturbacion de campo, suponiendo que esa matriz
# representa al litio


Chi =  24.1*1e-6 #(ppm) Susceptibilidad volumetrica

litio = matriz_3D * Chi
voxelsize = [0.01,0.01,0.01] #mm

delta = cFS.calculateFieldShift(litio, voxelsize)*1e6

#%% GRAFICOS

# N = [256,256,256]

# slice en z
z_slice = int(N[0]/2)
plt.figure(2)
plt.pcolormesh(matriz_3D[z_slice,:,:])
plt.xlabel('x')
plt.ylabel('y')

v = 5
plt.figure(20)
plt.pcolormesh(delta[z_slice,:,:], cmap='seismic', vmax=v, vmin=-v)
# contour son curvas de nivel, probalo:
#plt.contour(delta[z_slice,:,:], 50, cmap='seismic', vmax=v, vmin=-v)
plt.xlabel('x')
plt.ylabel('y')
plt.colorbar()

# slice en x
x_slice = int(N[2]/2)
plt.figure(3)
plt.pcolormesh(matriz_3D[:,:,x_slice])
plt.xlabel('y')
plt.ylabel('z')

v = 4
plt.figure(30)
plt.pcolormesh(delta[:,:,x_slice], cmap='seismic', vmax=v, vmin=-v)
# contour son curvas de nivel, probalo:
#plt.contour(delta[:,:,x_slice], 50, cmap='seismic', vmax=v, vmin=-v)
plt.xlabel('y')
plt.ylabel('z')
plt.colorbar()

#ver la Delta únicamente en función de z
plt.figure(40)
x = np.arange(0,256,1)
g = delta[:,int(N[1]/2),int(N[2]/2)]
h = delta[:,100,100]
ax = plt.subplot(111)
ax.plot(x,g,'-',linewidth=2, label=' x=y=128 ')
ax.plot(x,h,'-',linewidth=2, label=' x=y=100 ')
ax.set_ylim(-12, 6)
plt.xlabel('z [voxel]')
plt.ylabel('Delta [ppm]')
plt.title('Delta en función de z')
ax.legend()


#%%


#%%
  # al final del codigo, para que se vean los gradicos hay que poner
plt.show()
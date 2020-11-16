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

z_slice = int(N[0]/2)-int(nsz/2)
y_dim = np.arange(0,25.6,0.1)
x_dim = np.arange(0,25.6,0.1)


plt.figure(10)
plt.pcolormesh(x_dim,y_dim,matriz_3D[z_slice,:,:])
plt.xlabel('x [mm]')
plt.ylabel('y [mm]')


v = 5
plt.figure(11)

plt.pcolormesh(x_dim,y_dim,delta[z_slice,:,:], cmap='seismic', vmax=v, vmin=-v)
plt.xlabel('x [mm]')
plt.ylabel('y [mm]')
plt.colorbar()

# slice en x

x_slice = int(N[2]/2)
z_dim = np.arange(0,7.68,0.03)
y_dim = np.arange(0,25.6,0.1)

plt.figure(20)
plt.pcolormesh(y_dim,z_dim,matriz_3D[:,:,x_slice])
plt.xlabel('y [mm]')
plt.ylabel('z [mm]')


v = 4
plt.figure(21)
plt.pcolormesh(y_dim,z_dim,delta[:,:,x_slice], cmap='seismic', vmax=v, vmin=-v)
plt.xlabel('y [mm]')
plt.ylabel('z [mm]')
plt.colorbar()

#%%
#GRAFICOS PARA EL ANALISIS
#Delta en función de z para distintos radios R
#Lo hago dejando la posición en y fija y me paro en 4 valores distintos
#a lo largo de x, entonces tengo 4 pares (y,x). Esto es lo mismo que
#un análisis radial, dada la simetría de la muestra.


plt.figure(30)
x = np.arange(0,7.68,0.03)
f = delta[:,int(N[1]/2),int(N[2]/2)]
g = delta[:,int(N[1]/2),int(N[2]*79/128)]
h = delta[:,int(N[1]/2),int(N[2]*173/256)]
i = delta[:,int(N[1]/2),int(N[2]*47/64)]
ax = plt.subplot(111)
ax.plot(x,f,'-',linewidth=2, label=' R = 0 mm')
ax.plot(x,g,'-',linewidth=2, label=' R = 3 mm')
ax.plot(x,h,'-',linewidth=2, label=' R = 4.5 mm')
ax.plot(x,i,'-',linewidth=2, label=' R = 6 mm')
ax.set_ylim(-12, 7)
plt.xlabel('z [mm]')
plt.ylabel('Delta [ppm]')
plt.title('Delta(z)')
ax.legend()

#Veo unicamente la variación dentro de la celda para d_in
plt.figure(31)
x = np.arange(3.45,4.26,0.03)
z = np.arange(115,142,1)
f = delta[z,int(N[1]/2),int(N[2]/2)]
g = delta[z,int(N[1]/2),int(N[2]*79/128)]
h = delta[z,int(N[1]/2),int(N[2]*173/256)]
i = delta[z,int(N[1]/2),int(N[2]*47/64)]
ax = plt.subplot(111)
ax.plot(x,f,'-',linewidth=2, label=' R = 0 mm')
ax.plot(x,g,'-',linewidth=2, label=' R = 3 mm')
ax.plot(x,h,'-',linewidth=2, label=' R = 4.5 mm')
ax.plot(x,i,'-',linewidth=2, label=' R = 6 mm')
ax.set_ylim(-12, 7)
plt.xlabel('z [mm]')
plt.ylabel('Delta [ppm]')
plt.title('Delta(z) dentro de la celda')
ax.legend()

#Veo unicamente la variación fuera de la celda para d_out
plt.figure(32)
x = np.arange(4.17,7.68,0.03)
z = np.arange(139,256,1)
f = delta[z,int(N[1]/2),int(N[2]/2)]
g = delta[z,int(N[1]/2),int(N[2]*79/128)]
h = delta[z,int(N[1]/2),int(N[2]*173/256)]
i = delta[z,int(N[1]/2),int(N[2]*47/64)]
ax = plt.subplot(111)
ax.plot(x,f,'-',linewidth=2, label=' R = 0 mm')
ax.plot(x,g,'-',linewidth=2, label=' R = 3 mm')
ax.plot(x,h,'-',linewidth=2, label=' R = 4.5 mm')
ax.plot(x,i,'-',linewidth=2, label=' R = 6 mm')
ax.set_ylim(-12, 7)
plt.xlabel('z [mm]')
plt.ylabel('Delta [ppm]')
plt.title('Delta(z) fuera de la celda')
ax.legend()



#%%
  # al final del codigo, para que se vean los gradicos hay que poner
plt.show()
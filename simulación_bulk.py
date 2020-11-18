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

#Voy a construir un electrodo cilindrica, con las dimensiones h=0.71mm
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

# Defino los indices que conforman el electrodo 

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

#Veo unicamente la variación dentro del electrodo para Delta_in
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
plt.title('Delta(z) dentro del electrodo')
ax.legend()

#Veo unicamente la variación fuera del electrodo para Delta_out
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
plt.title('Delta(z) fuera del electrodo')
ax.legend()

#Gráficos para figuras esquemáticas
#Delta_in

plt.figure(33)
x = np.arange(3.45,4.26,0.03)
z = np.arange(115,142,1)
f = delta[z,int(N[1]/2),int(N[2]/2)]
ax = plt.subplot(111)
ax.plot(x,f,'-',linewidth=2, label=' R = 0 mm')
ax.set_ylim(-12, 7)
plt.xlabel('z [mm]')
plt.ylabel('Delta [ppm]')
plt.title('Delta(z) dentro del electrodo')
ax.legend()

#Delta_out
plt.figure(34)
x = np.arange(4.17,7.68,0.03)
z = np.arange(139,256,1)
f = delta[z,int(N[1]/2),int(N[2]/2)]
ax = plt.subplot(111)
ax.plot(x,f,'-',linewidth=2, label=' R = 0 mm')
ax.set_ylim(-12, 7)
plt.xlabel('z [mm]')
plt.ylabel('Delta [ppm]')
plt.title('Delta(z) fuera del electrodo')
ax.legend()

#Delta_centro

plt.figure(35)
x = np.arange(0,7.68,0.03)
f = delta[:,int(N[1]/2),int(N[2]/2)]
h = delta[:,int(N[1]/2),int(N[2]*173/256)]
ax = plt.subplot(111)
ax.plot(x,f,'-',linewidth=2, label=' R = 0 mm')
ax.plot(x,h,'-',linewidth=2, label=' R = 4.5 mm')
ax.set_ylim(-12, 7)
plt.xlabel('z [mm]')
plt.ylabel('Delta [ppm]')
plt.title('Delta(z)')
ax.legend()

#%%
#En esta sección determinamos los valores exactos para calcular Delta_in /
#Delta_out y Delta_centro
#NOTACIÓN: para clarificar los máximos y mínimos, iran con el subindice
# de su función seguido de los números 1 correspondiente a Delta_in,
# 2 a Delta_out y Delta_centro se calcula con estos valores

#Delta_in
#min
z = np.arange(115,142,1)
f = delta[z,int(N[1]/2),int(N[2]/2)]
g = delta[z,int(N[1]/2),int(N[2]*79/128)]
h = delta[z,int(N[1]/2),int(N[2]*173/256)]
i = delta[z,int(N[1]/2),int(N[2]*47/64)]

minf1 = np.min(f)
print(minf1,['minf1'])
ming1 = np.min(g)
print(ming1,['ming1'])
minh1 = np.min(h)
print(minh1,['minh1'])
mini1 = np.min(i)
print(mini1,['mini1'])
#max
z = np.arange(120,136,1)
f = delta[z,int(N[1]/2),int(N[2]/2)]
g = delta[z,int(N[1]/2),int(N[2]*79/128)]
h = delta[z,int(N[1]/2),int(N[2]*173/256)]
i = delta[z,int(N[1]/2),int(N[2]*47/64)]

maxf1 = np.max(f)
print(maxf1,['maxf1'])
maxg1 = np.max(g)
print(maxg1,['maxg1'])
maxh1 = np.max(h)
print(maxh1,['maxh1'])
maxi1 = np.max(i)
print(maxi1,['maxi1'])

#Diferencias
Delta_in_f = maxf1-minf1 
print(Delta_in_f,['Delta_in_f'])
Delta_in_g = maxg1-ming1
print(Delta_in_g,['Delta_in_g'])
Delta_in_h = maxh1-minh1
print(Delta_in_h,['Delta_in_h'])
Delta_in_i = maxi1-mini1
print(Delta_in_i,['Delta_in_i'])

#Delta_out
#min
z = np.arange(142,256,1)
f = delta[z,int(N[1]/2),int(N[2]/2)]
g = delta[z,int(N[1]/2),int(N[2]*79/128)]
h = delta[z,int(N[1]/2),int(N[2]*173/256)]
i = delta[z,int(N[1]/2),int(N[2]*47/64)]

minf2 = np.min(f)
print(minf2,['minf2'])
ming2 = np.min(g)
print(ming2,['ming2'])
minh2 = np.min(h)
print(minh2,['minh2'])
mini2 = np.min(i)
print(mini2,['mini2'])
#max
z = np.arange(139,256,1)
f = delta[z,int(N[1]/2),int(N[2]/2)]
g = delta[z,int(N[1]/2),int(N[2]*79/128)]
h = delta[z,int(N[1]/2),int(N[2]*173/256)]
i = delta[z,int(N[1]/2),int(N[2]*47/64)]

maxf2 = np.max(f)
print(maxf2,['maxf2'])
maxg2 = np.max(g)
print(maxg2,['maxg2'])
maxh2 = np.max(h)
print(maxh2,['maxh2'])
maxi2 = np.max(i)
print(maxi2,['maxi2'])

#Diferencias
Delta_out_f = maxf2-minf2 
print(Delta_out_f,['Delta_out_f'])
Delta_out_g = maxg2-ming2
print(Delta_out_g,['Delta_out_g'])
Delta_out_h = maxh2-minh2
print(Delta_out_h,['Delta_out_h'])
Delta_out_i = maxi2-mini2
print(Delta_out_i,['Delta_out_i'])



#%%
 
plt.show()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 07:22:30 2020

@author: muri
"""


import numpy as np 
import matplotlib.pyplot as plt 
from oct2py import Oct2Py
from mpl_toolkits.mplot3d import Axes3D 
import calculateFieldShift as cFS 
from SimulationVolume import *

#Voy a construir un electrodo cilindrica, con las dimensiones h=0.71mm
# y d=12mm. Además quiero FOVz~10*h y FOVxy~2*d


N = [512,256,256]
# Sólo por las dimensiones pedidas y el nro elegido de N, tomo este voxelsize
voxelSize = [0.015, 0.1, 0.1]# mm
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

# N = [512,256,256]

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
plt.clim(-4,4)
plt.colorbar()

# slice en x

x_slice = int(N[2]/2)
z_dim = np.arange(0,7.68,0.015)
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
plt.clim(-4,4)
plt.colorbar()

# Agrego un gráfico 3D del bulk
fig = plt.figure(60)
ax = fig.gca(projection='3d')
x = np.linspace(0,99,100)
ax.voxels(matriz_3D, facecolors='grey', edgecolor='k')

#%%
#GRAFICOS PARA VERIFICAR LA EXPORTACIÓN DE DATA DE PERFILES
#Delta en función de z para distintos radios R
#Lo hago dejando la posición en y fija y me paro en 4 valores distintos
#a lo largo de x, entonces tengo 4 pares (y,x). Esto es lo mismo que
#un análisis radial, dada la simetría de la muestra.


plt.figure(30)
x = np.arange(0,7.68,0.015)
f = delta[:,int(N[1]/2),int(N[2]/2)]
g = delta[:,int(N[1]/2),int(N[2]*79/128)]
h = delta[:,int(N[1]/2),int(N[2]*173/256)]
i = delta[:,int(N[1]/2),int(N[2]*186/256)]
ax = plt.subplot(111)
ax.plot(x,f,'-',linewidth=2, label=' R = 0 mm')
ax.plot(x,g,'-',linewidth=2, label=' R = 3 mm')
ax.plot(x,h,'-',linewidth=2, label=' R = 4,5 mm')
ax.plot(x,i,'-',linewidth=2, label=' R = 5,8 mm')
ax.set_ylim(-10, 10)
plt.xlabel('z [mm]')
plt.ylabel(r'$\delta$(z) [ppm]')
plt.title(' ')
ax.legend()


#%%
#Exportación de datos para el código de superposición
#PRUEBA

#R=0

#Grafico de la data a exportar
plt.figure(40)
x = np.arange(3.84,7.68,0.015)
z = np.arange(256,512,1)
f = delta[z,int(N[1]/2),int(N[2]/2)]
ax = plt.subplot(111)
ax.plot(x,f,'.',linewidth=2, label=' R = 0 mm')
ax.set_ylim(-12, 8)
plt.xlabel('z [mm]')
plt.ylabel('Delta [ppm]')
plt.title('Delta(z)')
ax.legend()

#Exportación de la data partida en z0

z0 = 280 #voxel

z_in = np.arange(256,int(z0),1)
z_out = np.arange(int(z0+1),512,1)
f_in = delta[z_in,int(N[1]/2),int(N[2]/2)]
f_out = delta[z_out,int(N[1]/2),int(N[2]/2)]

datos_in = np.array([z_in,f_in]).T
np.savetxt('perfil_radio3840.in',datos_in)

datos_out = np.array([z_out,f_out]).T
np.savetxt('perfil_radio3840.out',datos_out)


#hago una purueba para ver si los datos expordatos estan bien

#datos_cargados tiene la data de z está en la primer columna y la data de
# delta en la segunda columna

def loadtext_without_columns(filePath, skipcols=[], delimiter=" "):

    with open(filePath) as f:
 
        n_cols = len(f.readline().split(delimiter))

    #define a range from 0 to n_cols
    usecols = np.arange(0, n_cols)

    #remove the indices found in skipcols
    usecols = set(usecols) - set(skipcols)

    #sort the new indices in ascending order
    usecols = sorted(usecols)

    #load the file and retain indices found in usecols
    data = np.loadtxt('perfil_radio3840.in', usecols = usecols)

    return data
data_z = loadtext_without_columns("./perfil_radio3840.in",skipcols = [1], delimiter = " ")
print(data_z)
data_f = loadtext_without_columns("./perfil_radio3840.in",skipcols = [0], delimiter = " ")
print(data_f)

plt.figure(23)
ax = plt.subplot(111)
ax.plot(data_z,data_f,'o',linewidth=8, label=' datos in')
ax.plot(z,f,'.',linewidth=2, label=' datos originales')
ax.set_ylim(-12, 8)
plt.xlabel('z [voxel]')
plt.ylabel('Delta [ppm]')
plt.title('Delta')
ax.legend()

#%%
#AHORA SÍ EXPORTAMOS LOS DATOS CORRESPONDIENTES A LOS 4 RADIOS

#IMPORTANTE: El último voxel con información in es el 279 y el primer voxel
#con información out es el 280. Esto vale para todas las funciones de los 
#radios, se demuestran estos valores para f, pero se utilizan los mismos
# para g, h e i


#Grafico de la data a exportar
plt.figure(50)
x = np.arange(256,512,1)
z = np.arange(256,512,1)
f = delta[z,int(N[1]/2),int(N[2]/2)]
g = delta[z,int(N[1]/2),int(N[2]*79/128)]
h = delta[z,int(N[1]/2),int(N[2]*173/256)]
i = delta[z,int(N[1]/2),int(N[2]*186/256)]
ax = plt.subplot(111)
ax.plot(x,f,'.',linewidth=2, label=' R = 0 mm')
ax.plot(x,g,'.',linewidth=2, label=' R = 3 mm')
ax.plot(x,h,'.',linewidth=2, label=' R = 4.5 mm')
ax.plot(x,i,'.',linewidth=2, label=' R = 5.8 mm')
ax.set_ylim(-12, 10)
plt.xlabel('z [mm]')
plt.ylabel('Delta [ppm]')
plt.title('Delta(z)')
ax.legend()

#*****************************************************************************
#R=0mm corresponte realmente a 12800um

#Determinación de z0 para R=0

z = np.arange(3.84,7.68,0.015)
z_dom = np.arange(256,512,1)
f_det = delta[z_dom,int(N[1]/2),int(N[2]/2)]

min_f_det = np.min(f_det)
max_f_det = np.max(f_det)
z0_in = np.where(f_det == min_f_det)
z0_out = np.where(f_det == max_f_det)
print(z0_in)
print(z0_out)
z0_in = 4.185
print(z0_in)

#Exportación de la data partida en z0_0
 
z_in = np.arange(-0,3.84-z0_in,-0.015)
z_in_dom = np.arange(279,256,-1)
z_out = np.arange(0,3.48,0.015)
z_out_dom = np.arange(280,512,1) 
f_in = delta[z_in_dom,int(N[1]/2),int(N[2]/2)]
f_out = delta[z_out_dom,int(N[1]/2),int(N[2]/2)]

datos_in = np.array([z_in,f_in]).T
np.savetxt('perfil_radio000.in',datos_in)

datos_out = np.array([z_out,f_out]).T
np.savetxt('perfil_radio000.out',datos_out)

#*****************************************************************************
#R=3mm corresponde realmente a 15800um
#Exportación de la data partida en z0_3
 
z_in = np.arange(-0,3.84-z0_in,-0.015)
z_in_dom = np.arange(279,256,-1)
z_out = np.arange(0,3.48,0.015)
z_out_dom = np.arange(280,512,1) 
g_in = delta[z_in_dom,int(N[1]/2),int(N[2]*79/128)]
g_out = delta[z_out_dom,int(N[1]/2),int(N[2]*79/128)]

datos_in = np.array([z_in,g_in]).T
np.savetxt('perfil_radio300.in',datos_in)

datos_out = np.array([z_out,g_out]).T
np.savetxt('perfil_radio300.out',datos_out)

#*****************************************************************************
#R=4.5mm corresponde realmente a 16800um
#Exportación de la data partida en z0_4
 
z_in = np.arange(-0,3.84-z0_in,-0.015)
z_in_dom = np.arange(279,256,-1)
z_out = np.arange(0,3.48,0.015)
z_out_dom = np.arange(280,512,1) 
h_in = delta[z_in_dom,int(N[1]/2),int(N[2]*173/256)]
h_out = delta[z_out_dom,int(N[1]/2),int(N[2]*173/256)]

datos_in = np.array([z_in,h_in]).T
np.savetxt('perfil_radio450.in',datos_in)

datos_out = np.array([z_out,h_out]).T
np.savetxt('perfil_radio450.out',datos_out)

#*****************************************************************************
#R=5.8mm corresponde realmente a 18600um
#Exportación de la data partida en z0_3
 
z_in = np.arange(-0,3.84-z0_in,-0.015)
z_in_dom = np.arange(279,256,-1)
z_out = np.arange(0,3.48,0.015)
z_out_dom = np.arange(280,512,1) 
i_in = delta[z_in_dom,int(N[1]/2),int(N[2]*186/256)]
i_out = delta[z_out_dom,int(N[1]/2),int(N[2]*186/256)]

datos_in = np.array([z_in,i_in]).T
np.savetxt('perfil_radio580.in',datos_in)

datos_out = np.array([z_out,i_out]).T
np.savetxt('perfil_radio580.out',datos_out)
#%%
 
plt.show()
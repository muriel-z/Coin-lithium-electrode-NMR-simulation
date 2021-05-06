#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 07:22:30 2020

@author: muri
"""


import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.colors
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

z_slice = int(N[0]/2)-int(nsz/2)-1
y_dim = np.arange(-12.8,12.8,0.1)
x_dim = np.arange(-12.8,12.8,0.1)


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
z_dim = np.arange(-3.84,3.84,0.015)
y_dim = np.arange(-12.8,12.8,0.1)

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
# #%%
# # Agrego un gráfico 3D del bulk
# tmpvol =np.zeros((Nmz+5,Nmy,Nmx))
# tmpvol[1:-4,:,:] = matriz_3D
# tmpvol[0,:,:] = 1
# filename = './tmp.stl'
# with Oct2Py() as oc:
#   print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++")
#   print("Creando figura 3D. Esto puede demorar varios minutos...")
#   fv = oc.isosurface(tmpvol, 0.5) # Make patch w. faces "out"
#   oc.stlwrite(filename,fv)        # Save to binary .stl
# print("       Listo!") 
# print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++")

#%%
muestra = matriz_3D.T
fig = plt.figure(60)
ax = fig.gca(projection='3d')
ax.voxels(muestra, facecolors='grey', edgecolor='grey')


#%% 
#Grafico 3D con distintos colores y sombras
muestra = matriz_3D.T
def midpoints(x):
    sl = ()
    for i in range(x.ndim):
        x = (x[sl + np.index_exp[:-1]] + x[sl + np.index_exp[1:]]) / 2.0
        sl += np.index_exp[:]
    return x

# prepare some coordinates, and attach rgb values to each
r, g, b = np.indices((257, 257, 513)) / 512.0
rc = midpoints(r)
gc = midpoints(g)
bc = midpoints(b)
# combine the color components
colors = np.zeros(muestra.shape + (3,))
colors[..., 0] = rc
colors[..., 1] = gc
colors[..., 2] = bc

x_p,y_p,z_p = np.indices((257,257,513))

# and plot everything
ax = plt.figure(63).add_subplot(projection='3d')
ax.voxels(x_p,y_p,z_p, muestra,
          facecolors=colors,
          edgecolors=np.clip(2*colors - 0.5, 0, 1),  # brighter
          linewidth=0.5)
ax.set(xlabel='x [voxels]', ylabel='y [voxels]', zlabel='z [voxels]')
#%%
#Grafico para verificar la homogeneidad de delta(x), o sea perfil en x

plt.figure(70)
x_int = np.arange(69,188,1)
x = np.arange(-5.9,6.0,0.1)
f = delta[int(N[0]/2)+int(nsz/2)+1,int(N[1]/2),x_int]
ax = plt.subplot(111)
ax.plot(x,f,linestyle='-', marker='.',linewidth=2, label=r'$\delta$(x) ')
ax.set_xlim(-6,6)
plt.xlabel('x [mm]')
plt.ylabel(r'$\delta$(x) [ppm]')
plt.title(' ')
ax.legend()

#%%
#Grafico de comparación R=0 y tita de heavyside en las proximidades ocurrentes
plt.figure(71)
z = np.arange(273,290,1)
x = np.arange(0.255,0.51,0.015)
f = delta[z,int(N[1]/2),int(N[2]/2)]
i = delta[z,int(N[1]/2),int(N[2]*186/256)]
g = 15.97*np.heaviside(f,1)-8.60
ax = plt.subplot(111)
ax.plot(x,g,'-',linewidth=3,color='#c8d11f', label= 'Función escalón')
ax.plot(x,f,linestyle='-', marker='.', linewidth=2,color='#1f77b4', label=' Perfil R = 0 mm')
ax.plot(x,i,linestyle='-', marker='.', linewidth=2,color='#d62728', label=' Perfil R = 5,8 mm')
plt.xlabel('z [mm]')
plt.ylabel(r'$\delta$(z) [ppm]')
plt.title(' ')
ax.legend()
#%%
#Grafico para presentar eta_out 
plt.figure(72)
z = np.arange(279,290,1)
x = np.arange(0.345,0.50,0.015)
f = delta[z,int(N[1]/2),int(N[2]/2)]
ax = plt.subplot(111)
ax.plot(x,f,linestyle='-', marker='.',linewidth=2, label=' R = 0 mm')
ax.set_ylim(2,8)
plt.xlabel('z [mm]')
plt.ylabel(r'$\delta$(z) [ppm]')
plt.title(' ')
ax.legend()

#Grafico eta_in

plt.figure(73)
z = np.arange(273,281,1)
x = np.arange(0.255,0.375,0.015)
f = delta[z,int(N[1]/2),int(N[2]/2)]
ax = plt.subplot(111)
ax.plot(x,f,linestyle='-', marker='.',linewidth=2, label=' R = 0 mm')
ax.set_ylim(-9,-2)
plt.xlabel('z [mm]')
plt.ylabel(r'$\delta$(z) [ppm]')
plt.title(' ')
ax.legend()
#%%
#GRAFICOS PARA VERIFICAR LA EXPORTACIÓN DE DATA DE PERFILES
#Delta en función de z para distintos radios R
#Lo hago dejando la posición en y fija y me paro en 4 valores distintos
#a lo largo de x, entonces tengo 4 pares (y,x). Esto es lo mismo que
#un análisis radial, dada la simetría de la muestra.


plt.figure(30)
x = np.arange(-3.84,3.84,0.015)
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
#Agrego la seccion horrible para calcular los valores máximos y minimos de eta_in y eta_out

#En esta sección determinamos los valores exactos para calcular Delta_in /
#Delta_out y Delta_centro
#NOTACIÓN: para clarificar los máximos y mínimos, iran con el subindice
# de su función seguido de los números 1 correspondiente a eta_in,
# 2 a eta_out 

#eta_in
#min
z = np.arange(273,279,1)
f = delta[z,int(N[1]/2),int(N[2]/2)]
g = delta[z,int(N[1]/2),int(N[2]*79/128)]
h = delta[z,int(N[1]/2),int(N[2]*173/256)]
i = delta[z,int(N[1]/2),int(N[2]*186/256)]

minf1 = np.min(f)
print(minf1,['minf1'])
ming1 = np.min(g)
print(ming1,['ming1'])
minh1 = np.min(h)
print(minh1,['minh1'])
mini1 = np.min(i)
print(mini1,['mini1'])
#max
z = np.arange(273,279,1)
f = delta[z,int(N[1]/2),int(N[2]/2)]
g = delta[z,int(N[1]/2),int(N[2]*79/128)]
h = delta[z,int(N[1]/2),int(N[2]*173/256)]
i = delta[z,int(N[1]/2),int(N[2]*186/256)]

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

#eta_out
#min
z = np.arange(280,290,1)
f = delta[z,int(N[1]/2),int(N[2]/2)]
g = delta[z,int(N[1]/2),int(N[2]*79/128)]
h = delta[z,int(N[1]/2),int(N[2]*173/256)]
i = delta[z,int(N[1]/2),int(N[2]*186/256)]

minf2 = np.min(f)
print(minf2,['minf2'])
ming2 = np.min(g)
print(ming2,['ming2'])
minh2 = np.min(h)
print(minh2,['minh2'])
mini2 = np.min(i)
print(mini2,['mini2'])
#max
z = np.arange(280,290,1)
f = delta[z,int(N[1]/2),int(N[2]/2)]
g = delta[z,int(N[1]/2),int(N[2]*79/128)]
h = delta[z,int(N[1]/2),int(N[2]*173/256)]
i = delta[z,int(N[1]/2),int(N[2]*186/256)]

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
 



plt.show()
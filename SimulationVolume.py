#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 12:45:24 2020

@author: santi
"""

import numpy as np

def SimulationVolume(N=None, voxelSize=None, FOV=None):
  """
  Esta funcion determina el volumen de la simulacion.
  es decir determina N, voxelSize o FOV, alguno de ellos dados los otros 
  dos restantes.
  
  N, voxelSize y FOV deben ser array_like, de 3 elementos con la convencion
  [z,y,x]
  
  ejemplo: dado N y FOV (en mmm) determino voxelSize
    
  N, voxelSize, FOV = SimulationVolume(N=[256,64,64], FOV=[1.5,10,10])
  
  Notar que el resultado son los tres parametros, ya que de poso los transforma
  en np.arrays!
  """ 
  if FOV is None:
    N = np.array(N)
    voxelSize = np.array(voxelSize)  
    FOV = N*voxelSize
  elif voxelSize is None:  
    N = np.array(N)
    FOV = np.array(FOV)
    voxelSize = FOV/N.astype(float)
  elif N is None:  
    FOV = np.array(FOV)
    voxelSize = np.array(voxelSize)  
    N = FOV/voxelSize
    
  N = N.astype(int)
  #mensaje en pantalla:
  sN = str(N[0])+'*'+str(N[1])+'*'+str(N[2])
  sFOV = str(FOV[0])+'mm*'+str(FOV[1])+'mm*'+str(FOV[2])+'mm'
  sVS = str(voxelSize[0])+'mm*'+str(voxelSize[1])+'mm*'+str(voxelSize[2])+'mm'
  print('Volumen de simulacion (convencion z,y,x):')
  print('N: '+ sN + ' | voxelSize:  ' + sVS + ' | FOV:  '+ sFOV)
  return N, voxelSize, FOV  
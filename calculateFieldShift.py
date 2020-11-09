"""
Created on Tue Dec 17 12:15:20 2019

@author: santi

Voy a "tradicir" el codigo de calculateFieldShift.m a python
"""
import numpy as np
from scipy import fftpack
import time

# Fourier Based Field-shift Calculation for MRI
# - k-space Filtering (KF) implementation
# - implemented with 'virtual' zero-padding
# Written 22-06-2012 by Job Bouwman (jgbouwman@hotmail.com)
# Update 14-06-2013: to cope with anisotropic voxel size
# Update 26-01-2014: implemented a numerical trick to further increase the 
#           calculation speed with ~30% . As all distributions are real
#           in the spatial domain, the imaginary channel of the convolution
#           can be used for a second convolution. This channel is exploited 
#           to calculate the aliasing distribution.
# Update 20-10-2014: the use of nested functions allows larger matrices
#           again 13% faster, by using more efficient nested functions 
# method description (including scientific references)
# This is an efficient implementation of the Fourier-based procedure to
# calculate the susceptibility-induced field-shift if an object 
# with known susceptibility distribution (dChi_3D) is placed within 
# a strong homogeneous magnetic field (B0).   
# This method is based on the publication of: 
# Salomir et al, "A fast calculation method for magnetic field 
# inhomogeneity due to an arbitrary distribution of bulk susceptibility"
# Concepts in Magnetic Resonance Part B: Magnetic Resonance Engineering, 
# Volume 19B, Issue 1, pages 26?34, 2003
# However, the original method needs (inefficient zero-padding)
# Cheng YCN, Neelavalli J, Haacke EM. "Limitations of calculating field
# distributions and magnetic susceptibilities in MRI using a Fourier
# based method". Phys Med Biol 2009;54:1169?1189
# In the proposed method here, aliasing is removed by estimating 
# this artifact using an additional convolution in a lower resolution, 
# which is feasible due to the smooth decay of the dipole function. 
# The article describing the principle of low resolution aliasing 
# subtraction:
# Job G. Bouwman, Chris J.G. Bakker. "Alias subtraction more efficient than
# conventional zero-padding in the Fourier-based calculation of the 
# susceptibility induced perturbation of the magnetic field in MR" 
# 2012, Magnetic Resonance in Medicine 68:621630 (2012). 
# Please refer to this article when using this script.
# How to use this file:
# dChi_3D is your 3D susceptibility distribution
#   * the third dimension of dChi_3D is the one parallel to B0.
#   * the values of your dChi_3D are assumed to be relative to the 
#     dChi of the embedding medium. In practice this will either be: 
#     - relative to the dChi of air (0.36 ppm) if you want to calculate the 
#       field for example in a complete head model
#     - relative to the dChi of tissue (-9.05 ppm) if you want to calculate
#       the local effect of for example a small implant embedded in tissue.
# voxelSize (either in mm, centimeter, or inches)
#    * we use the matlabconvention: [dy, dx, dz]
# Example:
#   B0 = 3; % tesla;
#   voxelSize = [1 1 2]; % mm
#   gyromagneticRatio = 2*pi*42.58e6
#   dB = B0*calculateFieldShift(dChi_3D, voxelSize);
#   dOhmega = dB*gyromagneticRatio;

def calculateFieldShift(dChi_3D, voxelSize):
    """
    voxelSize tiene que ser un array de tres elementos 
    """
    #     Checking the input parameters:
    #    if or(nargin < 1, nargin > 2)
    #        error('Please supply two arguments (dChi_3D and voxelSize).');
    #    end
    #    dChi_3D = varargin{1};
    #    if or(length(size(dChi_3D)) ~= 3,  min(size(dChi_3D)) < 2)
    #        error('Input dChi_3D must be of dimension 3.');
    #    end
    #    
    #    if nargin == 2
    #        voxelSize = squeeze(varargin{2});
    #        if or(size(voxelSize,2)~=3, size(voxelSize,1)~=1)
    #            error('voxelSize must be a vector of length 3.');
    #        end
    #    else
    #        display(['In calculateFieldShift.m the voxelSize was not '...
    #                 'specified. The script assumed it to be isotropic:']);
    #        voxelSize = [1,1,1]
    #    end        
    # making sure that the size of the distribution in each dimension 
    # is a multiple of four:
    #    [NyInput, NxInput, NzInput] = np.shape(dChi_3D);
    #    
    #    if NxInput%4 + NyInput%4 + NzInput%4 > 0
    #        conventionalPad = mod(4-[mod(NyInput,4),mod(NxInput,4),mod(NzInput,4)], 4);
    #        dChi_3D = padarray(dChi_3D,conventionalPad,'post');
    #    end
    
    
    NzInput, NyInput, NxInput = np.shape(dChi_3D)
    FOV = voxelSize*np.array([NzInput, NyInput, NxInput])
        
    # # The Fourier-based convolution:
    # 1) The forward Fourier Transform (FT):
    # NB: Instead of doing an extra convolution to calculate the 
    # effect of aliasing (as proposed in the article), here we use 
    # the 'imaginary channel' to do this. 
    # This is done by downscaling the artificial environment, and
    # adding this as an imaginary component to the original 
    # (unpadded) susceptibility distribution. Together this is 
    # referred to as the "DUAL" distribution:
    
    ZcX = Zc(D(dChi_3D))
    FT_dChi_3D_DUAL = fftpack.fftn(dChi_3D + 1j * ZcX) 
    del dChi_3D  # en matlab, borraba la variable
            
    # 2) multiplication with the dipole function: 
    FT_dField_3D = FT_dChi_3D_DUAL*KF_kernel(FOV,np.shape(FT_dChi_3D_DUAL))
    del FT_dChi_3D_DUAL
    # 3) inverse FT to spatial domain:
    dField_3D_DUAL = fftpack.ifftn(FT_dField_3D)
    del FT_dField_3D
    # 4) Subtracting the upscaled, cropped center of the aliasing:
    dField_3D = np.real(dField_3D_DUAL) - UC(np.imag(dField_3D_DUAL))
    del dField_3D_DUAL
  
    # returning the field shift in the same size:
    if NxInput%4 + NyInput%4 + NzInput%4 > 0:
        nx = np.arange(0,NxInput+1)
        ny = np.arange(0,NyInput+1)
        nz = np.arange(0,NzInput+1)
        dField_3D = dField_3D[nz,ny,nx]
    
    return dField_3D

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
    
def KF_kernel(FOV,N):
    """
    FOV : field of view in x, y, and z directions
    N   : no of samples in kx, ky, kz
    """
    
    x = np.arange(int(-N[2]/2),int(N[2]/2))/FOV[2]
    y = np.arange(int(-N[1]/2),int(N[1]/2))/FOV[1]
    z = np.arange(int(-N[0]/2),int(N[0]/2))/FOV[0]
    
    kx_squared = fftpack.ifftshift(x)**2;
    ky_squared = fftpack.ifftshift(y)**2;
    kz_squared = fftpack.ifftshift(z)**2;
    
    [kz2_3D,ky2_3D,kx2_3D] = np.meshgrid(kz_squared,ky_squared,kx_squared, indexing='ij')
    
    
    k2 = kx2_3D + ky2_3D + kz2_3D
    
    # no quiero que me tire advertencias la division
    old_settings = np.seterr(divide='ignore', invalid='ignore')
    # np.divide(<whatever>, 0.0) = inf
    kernel = 1/3 - np.divide(kz2_3D,k2)
    np.seterr(**old_settings)
    
    kernel[0,0,0] = 0
    
    return kernel
#------------------------------------------------------------------------------
def D(X):
    """
    downscales a 3D distribution by factor two, by collapsing each cube  
    of eight voxels into one single voxel (average value) 
    """
    DX =      X[ ::2, ::2, ::2]    
    DX = DX + X[1::2, ::2, ::2]    
    DX = DX + X[ ::2, ::2,1::2]    
    DX = DX + X[1::2, ::2,1::2]    
    DX = DX + X[ ::2,1::2, ::2]    
    DX = DX + X[1::2,1::2, ::2]
    DX = DX + X[ ::2,1::2,1::2]    
    DX = DX + X[1::2,1::2,1::2]
    DX = DX/8.
    
    return DX
#------------------------------------------------------------------------------
    
def Zc(X):
    """
    This operator is complementary to zero-padding: the result is 
    a void space embedded in an environment of copies of X. 
    The resulting 3D distribution is 2 times as large 
    as the original in each dimension
    """
    dim = np.shape(X)
    # numpy follows te convention zyx
    # So.. x is the 3rd dimension: dim[2]
    #      y is the 2nd dimension: dim[1]
    #      z is the 1st dimension: dim[0]
    Nx = 2*dim[2]
    Ny = 2*dim[1]
    Nz = 2*dim[0]
    mid_x = (int((Nx/4.)), int((Nx*3./4.))) # the slicing will be between
    mid_y = (int((Ny/4.)), int((Ny*3./4.))) # mid[0]:mid[1], which begins in 
    mid_z = (int((Nz/4.)), int((Nz*3./4.))) # mid[0] and ends in mid[1]-1

    #ZcX = repmat(circshift(X, [Ny/4,       Nx/4,       Nz/4]), [2 2 2]);  asi es en matlab
    ZcX = np.tile(np.roll  (X, [int(Nz/4.), int(Ny/4.), int(Nx/4.)], axis=[0,1,2]), (2,2,2))
    ZcX[mid_z[0]:mid_z[1], mid_y[0]:mid_y[1], mid_x[0]:mid_x[1]] = 0
    
    return ZcX

#------------------------------------------------------------------------------
def UC2(domain):
    """
    UC = Upscaling + Cropping:
    both crops as upscales (factor 2) the central part of a larger 
    distribution. Therefore the size of input ('domain') and output
    ('center') are the same.
    
    ESTE FUE EL PRIMER INTENTO
    """
    t = time.time()
    dim = np.shape(domain)        
    Nx = dim[2]
    Ny = dim[1]
    Nz = dim[0]
    

    Mx = np.round(np.arange(1,Nx+1)/2 + 1/4 + Nx/4).astype(int) - 1
    My = np.round(np.arange(1,Ny+1)/2 + 1/4 + Ny/4).astype(int) - 1
    Mz = np.round(np.arange(1,Nz+1)/2 + 1/4 + Nz/4).astype(int) - 1
  
    Sx = mx + (-1)**np.arange(1,Nx+1)
    Sy = my + (-1)**np.arange(1,Ny+1)
    Sz = mz + (-1)**np.arange(1,Nz+1)

    center = np.zeros_like(domain)
    center =          27./64.*domain[np.ix_(Mz,My,Mx)]
    center = center +  9./64.*domain[np.ix_(Sz,My,Mx)]
    center = center +  9./64.*domain[np.ix_(Mz,My,Sx)]
    center = center +  3./64.*domain[np.ix_(Sz,My,Sx)]
    center = center +  9./64.*domain[np.ix_(Mz,Sy,Mx)]
    center = center +  3./64.*domain[np.ix_(Sz,Sy,Mx)]
    center = center +  3./64.*domain[np.ix_(Mz,Sy,Sx)]
    center = center +  1./64.*domain[np.ix_(Sz,Sy,Sx)]
          
    elapsed = time.time() - t
    print('tiempo en UC:  ', elapsed)          
    return center
  
#------------------------------------------------------------------------------  
def UC(domain):
    """
    UC = Upscaling + Cropping:
    both crops as upscales (factor 2) the central part of a larger 
    distribution. Therefore the size of input ('domain') and output
    ('center') are the same.
    
    ESTE INTENTO USO np.ix_ en lugar del meschgrid, CON ESTO LOGRE BAJAR EL TIEMPO!
    IMPORTANTE: iniciar con el nucleo vacio
    """
    t = time.time()
    dim = np.shape(domain)        
    Nx = dim[2]
    Ny = dim[1]
    Nz = dim[0]
    

    mx = np.round(np.arange(1,Nx+1)/2 + 1/4 + Nx/4).astype(int) - 1
    my = np.round(np.arange(1,Ny+1)/2 + 1/4 + Ny/4).astype(int) - 1
    mz = np.round(np.arange(1,Nz+1)/2 + 1/4 + Nz/4).astype(int) - 1
  
    sx = mx + (-1)**np.arange(1,Nx+1)
    sy = my + (-1)**np.arange(1,Ny+1)
    sz = mz + (-1)**np.arange(1,Nz+1)
     
    
    Mz, My, Mx = np.meshgrid(mz,my,mx, indexing='ij')
    Sz, Sy, Sx = np.meshgrid(sz,sy,sx, indexing='ij')
    
    center = np.zeros_like(domain)
    center =          27./64.*domain[Mz,My,Mx]
    center = center +  9./64.*domain[Sz,My,Mx]
    center = center +  9./64.*domain[Mz,My,Sx]
    center = center +  3./64.*domain[Sz,My,Sx]
    center = center +  9./64.*domain[Mz,Sy,Mx]
    center = center +  3./64.*domain[Sz,Sy,Mx]
    center = center +  3./64.*domain[Mz,Sy,Sx]
    center = center +  1./64.*domain[Sz,Sy,Sx]
          
    elapsed = time.time() - t
    print('tiempo en UC:  ', elapsed)          
    return center
  
#------------------------------------------------------------------------------  
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------



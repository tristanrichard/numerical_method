# -*- coding: utf-8 -*-
"""
Created on Tue March 15 11:33:02 2025

@author: Tristan
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy as sc

def gaussian(x, sigma, U,kp):
    return U *np.cos(kp*x)* np.exp(-(x**2) / (sigma**2))

if __name__ == '__main__':
    ######Parameter####
    sigma=1
    U=1
    h=sigma/8
    L=16*sigma
    N=int(L/h)
    kp= (24*np.pi)/L

    x=np.linspace(-L/2,L/2,N, endpoint=False)

    ##### Fourrier coefficients and wave number#####
    F=sc.fft.fft(gaussian(x,sigma,U,kp),norm='forward') ##normalize with N
    F=sc.fft.fftshift(F)
    k = 2*np.pi*sc.fft.fftfreq(N, d=h)  
    k=sc.fft.fftshift(k)


    
    
    j=np.linspace(-N/2,N/2,N, endpoint=False)
    #### Figure ####
    plt.figure(figsize=(8, 6))
    plt.plot(j, np.log(np.abs(F)), label="Discrete Fourier transform", marker='o', linestyle='None',color="blue")
    plt.xlabel("Mode index j",fontsize=15)
    plt.ylabel(r"$\log(|\hat{F}_j|)$",fontsize=15)
    plt.tick_params(axis='both', labelsize=15)
    plt.ylim([-50,0])
    #plt.legend(fontsize=12,loc="upper left")
    plt.grid()
    plt.savefig('../plot/fft_wavepacket.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    
    k = (2*np.pi*j)/L
    E2 = np.cos(k*h)
    E4 = (((4/3)*np.cos(k*h))-((1/3)*(1-2*(np.sin(k*h)**2))))
    I4 = (((3/2)*np.cos(k*h))+(3/4))/((1+((1/2)*np.cos(k*h)))**2)
    ED = (-(1/3)*np.exp(-1j*2*k*h)) + (np.exp(-1j*k*h)) +((1/3)*np.exp(1j*k*h))
    plt.figure(figsize=(8, 6))
    plt.plot(j,E2, label="E2 scheme", linestyle='-',color="brown")
    plt.plot(j,E4, label="E4 scheme", linestyle='-',color='darkgreen')
    plt.plot(j,I4, label="I4 scheme", linestyle='-',color='darkorange')
    plt.plot(j,np.real(ED), label="ED scheme", linestyle='-',color='darkblue')
    plt.axhline(y=1,color='black',label='Perfect case',linestyle='--')


    plt.xlabel("Mode index j",fontsize=15)
    plt.ylabel(r"$\frac{c^*_g}{c}$",fontsize=25)
    plt.tick_params(axis='both', labelsize=15)
    plt.legend(fontsize=12,loc="lower center")
    plt.grid()
    plt.savefig('../plot/cg_wavepacket.png', dpi=300, bbox_inches='tight')
    plt.show()
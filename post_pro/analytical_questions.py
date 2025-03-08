# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 11:33:02 2025

@author: Tristan
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy as sc

def gaussian(x, sigma, U):
    return U * np.exp(-(x**2) / (sigma**2))

def k_star_h(kh):
    term1 = np.exp(-1j * kh * 2 )
    term2 = np.exp(-1j * kh )
    term3 = np.exp(1j * kh)
    
    k_star_h = (term1 - (6 * term2) + 3 + (2 * term3)) / (6j)
    return k_star_h

def rho(coef,CFL):
    x=coef*CFL
    return 1 + x + (0.5*(x**2)) + ((1/6)*(x**3)) + ((1/24)*(x**4))

def find_max_CFL(coef,typ):
    CFL=10
    step=0.001
    if typ=='max':
        while CFL>0:
            if np.all(np.abs(rho(coef,CFL)) <=1):
                return CFL
            else:
                CFL = CFL-step
                
    elif typ=='opt':
        while CFL>0:
            if np.all((0.95 < np.abs(rho(coef, CFL))) & (np.abs(rho(coef, CFL)) <= 1)):
                return CFL
            else:
                CFL = CFL-step
            
    return None       


if __name__ == '__main__':
    ######Parameter####
    sigma=1
    U=1
    domain_size="big" #or 'small'
    
    ##### Grid definition####
    h=sigma/8
    if(domain_size=='big'):
        L=16*sigma
        N=int(L/h)
    elif(domain_size=='small'):
        L=4*sigma
        N=int(L/h)
    x=np.linspace(-L/2,L/2,N, endpoint=False)
    print("######### grid configuration #########")
    print(f"Number of grid points:{N}")
    print(f"Grid Length : {L}")
    print(f'grid size: {h}')
    print(f'grid points: {x}')
    print("######################################\n")

    
    ##### Fourrier coefficients and wave number#####
    F=sc.fft.fft(gaussian(x,sigma,U),norm='forward') ##normalize with N
    F=sc.fft.fftshift(F)
    k = 2*np.pi*sc.fft.fftfreq(N, d=h)  
    k=sc.fft.fftshift(k)


    #### Real Fourrier ####
    F_real = U * np.sqrt(np.pi) * sigma * (1 / L) * np.exp(- ((sigma**2) * (k**2)) / 4)

    
    
    j=np.linspace(-N/2,N/2,N, endpoint=False)
    #### Figure ####
    plt.figure(figsize=(8, 6))
    plt.plot(j, np.log(np.abs(F)), label="Discrete Fourier transform", marker='o', linestyle='None',color="blue")
    plt.plot(j, np.log(np.abs(F_real)), label="Analytical Fourier transform", linestyle="dashed",color="black")
    plt.xlabel("Mode index j",fontsize=15)
    plt.ylabel(r"$\log(|\hat{F}_j|)$",fontsize=15)
    plt.tick_params(axis='both', labelsize=15)
    plt.ylim([-50,0])
    plt.legend(fontsize=12,loc="upper left")
    plt.grid()
    plt.savefig(f'plot/fft_{domain_size}.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    
    kh=k*h
    kstar=k_star_h(kh)##kstar_h in fact
    
    real_kstar = np.real(kstar)
    imag_kstar = np.imag(kstar)
    
    # Plotting the results
    plt.figure(figsize=(10, 8))
    plt.plot(kh/np.pi, real_kstar/np.pi, label='ED scheme', color='b', linestyle='-', linewidth=2)
    plt.plot(kh/np.pi,kh/np.pi,color='black',label='Perfect case')
    #plt.title(r'Real of $k^*h$ as a Function of $kh$', fontsize=14)
    plt.xlabel(r'$\frac{kh}{\pi}$', fontsize=25)
    plt.ylabel(r'$\frac{\Re(k^*)h}{\pi}$', fontsize=25)
    plt.legend(fontsize=12,loc='best')
    plt.grid(True)
    plt.xlim([0,1])
    plt.ylim([0,1])
    plt.tick_params(axis='both', labelsize=15)
    plt.savefig(f'plot/dispersion_error.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Plotting the results
    plt.figure(figsize=(10, 8))
    plt.plot(kh/np.pi, imag_kstar/np.pi, label='ED scheme', color='blue', linestyle='-', linewidth=2)
    plt.plot(kh/np.pi,kh*0,color='black',label='Perfect case')
    #plt.title(r'Imaginary of $k^*h$ as a Function of $kh$', fontsize=14)
    plt.xlabel(r'$\frac{kh}{\pi}$', fontsize=25)
    plt.ylabel(r'$\frac{\Im(k^*)h}{\pi}$', fontsize=25)
    plt.legend(fontsize=12,loc='best')
    plt.grid(True)
    plt.xlim([0,1])
    plt.ylim([-0.5,0.1])
    plt.tick_params(axis='both', labelsize=15)
    plt.savefig(f'plot/dissipation_error.png', dpi=300, bbox_inches='tight')
    plt.show()
        
    
    
    ##### CFL #####
    print('######### CFL #######') 

    CFL_ED=find_max_CFL( -1j * kstar,'max')
    print(f' CFL of the ED scheme is {CFL_ED}')
    
    ###E2###
    kstar_E2 = np.sin(kh)
    CFL_E2=find_max_CFL( -1j * kstar_E2,'max')
    print(f' CFL of the E2 scheme is {CFL_E2}')
    
    ###E4 ###
    kstar_E4 = np.sin(kh)*(1+((1/3)*(1-np.cos(kh))))
    CFL_E4=find_max_CFL( -1j * kstar_E4,'max')
    print(f' CFL of the E4 scheme is {CFL_E4}')
    
    ###I4####
    kstar_I4 = ((3/2)*np.sin(kh))/(1+((1/2)*np.cos(kh)))
    CFL_I4=find_max_CFL( -1j * kstar_I4,'max')
    print(f' CFL of the I4 scheme is {CFL_I4}')
    
    CFL_ED_opt=find_max_CFL( -1j * kstar,'opt')
    print(f' CFL opt of the ED scheme is {CFL_ED_opt}')
    
    ###E2###
    kstar_E2 = np.sin(kh)
    CFL_E2_opt=find_max_CFL( -1j * kstar_E2,'opt')
    print(f' CFL opt of the E2 scheme is {CFL_E2_opt}')
    
    ###E4 ###
    kstar_E4 = np.sin(kh)*(1+((1/3)*(1-np.cos(kh))))
    CFL_E4_opt=find_max_CFL( -1j * kstar_E4,'opt')
    print(f' CFL opt of the E4 scheme is {CFL_E4_opt}')
    
    ###I4####
    kstar_I4 = ((3/2)*np.sin(kh))/(1+((1/2)*np.cos(kh)))
    CFL_I4_opt=find_max_CFL( -1j * kstar_I4,'opt')
    print(f' CFL opt of the I4 scheme is {CFL_I4_opt}')
    



    
    
    
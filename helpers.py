import numpy as np

DENS = 0
MOMX = 1
MOMY = 2
ENER = 3
VELX = 1
VELY = 2
PRES = 3
FWD = 1
BWD = -1

def stabilization_roe(u_l, up_l, u_r, up_r):
    # Calculate left/right conservative variables
    #-----------------------------------------------------------------------------------------------------------------------------------
    # Calculate left/right enthalpy
    H_r=(u_r[:,ENER] + up_r[:,PRES])/u_r[:,DENS]
    H_l=(u_l[:,ENER] + up_l[:,PRES])/u_l[:,DENS]
    #-----------------------------------------------------------------------------------------------------------------------------------
    # Calculate sqrt(rho)
    sqrt_rho_r = np.sqrt(u_r[:,DENS])
    sqrt_rho_l = np.sqrt(u_l[:,DENS])
    sum_sqrt_rho_q = 1. / (sqrt_rho_l + sqrt_rho_r)
    #-----------------------------------------------------------------------------------------------------------------------------------
    # Calculate Roe mean values
    v1_bar=(sqrt_rho_r*up_r[:,VELX]+sqrt_rho_l*up_l[:,VELX]) * sum_sqrt_rho_q
    v2_bar=(sqrt_rho_r*up_r[:,VELY]+sqrt_rho_l*up_l[:,VELY]) * sum_sqrt_rho_q
    H_bar=(sqrt_rho_r* H_r+sqrt_rho_l* H_l) * sum_sqrt_rho_q
    c_bar=np.sqrt(0.4*(H_bar-0.5*(v1_bar*v1_bar + v2_bar*v2_bar)))
    #-----------------------------------------------------------------------------------------------------------------------------------
    # Calculate mean eigenvalues
    a1=v1_bar-c_bar
    a2=v1_bar
    a3=v1_bar
    a4=v1_bar+c_bar
    #-----------------------------------------------------------------------------------------------------------------------------------
    # Calculate mean eigenvectors
    r1 = np.zeros_like(u_l)
    r2 = np.zeros_like(u_l)
    r3 = np.zeros_like(u_l)
    r4 = np.zeros_like(u_l)

    r1[:,0]=1.
    r1[:,1]=a1
    r1[:,2]=v2_bar
    r1[:,3]=H_bar-v1_bar*c_bar

    r2[:,0]=1.
    r2[:,1]=v1_bar
    r2[:,2]=v2_bar
    r2[:,3]=0.5*(v1_bar*v1_bar + v2_bar*v2_bar)

    r3[:,0]=0.
    r3[:,1]=0.
    r3[:,2]=1.
    r3[:,3]=v2_bar

    r4[:,0]=1.
    r4[:,1]=a4
    r4[:,2]=v2_bar
    r4[:,3]=H_bar+v1_bar*c_bar
    #-----------------------------------------------------------------------------------------------------------------------------------
    # Calculate differences
    delta_rho = u_r[:,DENS] - u_l[:,DENS]
    delta_m1  = u_r[:,MOMX] - u_l[:,MOMX]
    delta_m2  = u_r[:,MOMY] - u_l[:,MOMY]
    delta_e   = u_r[:,ENER] - u_l[:,ENER]
    delta_eq  = delta_e-(delta_m2-v2_bar*delta_rho)*v2_bar
    #-----------------------------------------------------------------------------------------------------------------------------------
    # Calculate wave strengths gamma1,...,gamma4
    c_quer_q = 1./c_bar
    gamma2   = -0.4*c_quer_q*c_quer_q * (delta_rho*(v1_bar*v1_bar-H_bar)+ delta_eq - delta_m1*v1_bar )
    gamma1   = -0.5*c_quer_q*(delta_m1-delta_rho*(v1_bar+c_bar))-0.5*gamma2
    gamma4   = delta_rho - gamma1 - gamma2
    gamma3   = delta_m2 - v2_bar*delta_rho
    return   gamma1[...,np.newaxis]*np.abs(a1[...,np.newaxis])*r1 \
           + gamma2[...,np.newaxis]*np.abs(a2[...,np.newaxis])*r2 \
           + gamma3[...,np.newaxis]*np.abs(a3[...,np.newaxis])*r3 \
           + gamma4[...,np.newaxis]*np.abs(a4[...,np.newaxis])*r4
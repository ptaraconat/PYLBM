import matplotlib.pyplot as plt 
import sys as sys  
import numpy as np 
#LBM numerics

def GetLatice(D=2,Q=9):
    #ex = np.array([0,1,0,-1,0,1,-1,-1,1])
    #ey = np.array([0,0,1,0,-1,1,1,-1,-1])
    #w = np.array([4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36])
    ex = np.array([0,1,1,0,-1,-1,-1,0,1])
    ey = np.array([0,0,1,1,1,0,-1,-1,-1])
    w =np.array([4/9,1/9,1/36,1/9,1/36,1/9,1/36,1/9,1/36])
    cs = (1/np.sqrt(3)) 
    return ex, ey, w, cs

def CalcFeq(rho,ux,uy,ex,ey,w,cs,Nx,Ny,Q,thau,srct =[0,0]):
    ## Input 
    # rho[i,j] :::
    # ux[i,j] :::
    # uy[i,j] :::
    # ex[q] :::
    # ey[q] :::
    # w[q] :::
    # cs ::: speed of sound 
    # Nx :::
    # Ny :::
    # Q :::
    ## Output 
    # feq[i,j,q] :::
    #Init variables
    feq = np.zeros((Nx,Ny,Q))
    cs_squared = cs * cs 
    #print('Computing Feq')
    for i in range(Nx): 
        for j in range(Ny):
            ux[i,j] = ux[i,j] + srct[0]#(thau*srct[0]/rho[i,j])
            uy[i,j] = uy[i,j] + srct[1]#(thau*srct[1]/rho[i,j])
            Usq = ((ux[i,j]*ux[i,j])+(uy[i,j]*uy[i,j]))/cs_squared
            for q in range(Q):
                Ue = ((ex[q]*ux[i,j]) + (ey[q]*uy[i,j]))/cs_squared
                
                feq[i,j,q] = w[q]*rho[i,j]*(1 + Ue + 0.5*Ue*Ue -0.5*Usq )   
    return feq

def CalcMacro(f,ex,ey,Nx,Ny,Q):
    ## Input
    # f[i,j,q] :::
    # ex[q] :::
    # ey[q] :::
    # Nx :::
    # Ny :::
    # Q :::
    ## Output
    # rho[i,j] :::
    # ux[i,j] :::
    # uy[i,j] :::
    # Init variable 
    rho =  np.zeros((Nx,Ny))
    ux = np.zeros((Nx,Ny))
    uy = np.zeros((Nx,Ny))
    #print('Computing macroscopic quantities')
    # Calculate rho 
    for i in range(Nx): 
        for j in range(Ny):
            rho[i,j] = 0
            for q in range(Q):
                rho[i,j] = rho[i,j] + f[i,j,q]
    # Could be optimized by summing f along its third dimension
    # Calculate u
    for i in range(Nx): 
        for j in range(Ny):
            ux[i,j] = 0
            uy[i,j] = 0
            for q in range(Q):
                ux[i,j] = ux[i,j] + ex[q]*f[i,j,q]
                uy[i,j] = uy[i,j] + ey[q]*f[i,j,q]
            ux[i,j] = ux[i,j]/rho[i,j]
            uy[i,j] = uy[i,j]/rho[i,j]
            
    return ux,uy,rho

def Collision(f,feq,thau,Nx,Ny,Q,Dt=1,Bnd = None,with_boundaries = False):
    #Input 
    # f[i,j,q] :::
    # feq[i,j,q] :::
    # Dt :::
    # thau :::
    # Nx :::
    # Ny :::
    # Q :::
    #Output 
    # fstar[i,j,q]
    #Init variables
    if not with_boundaries :
        print('No boudary condition provided')
        Bnd = np.zeros((Nx,Ny))
    fstar = np.zeros((Nx,Ny,Q))
    #print('Collision step')
    for i in range(Nx): 
        for j in range(Ny):
            if (Bnd[i,j] == 0) or (Bnd[i,j] == 1):
                for q in range(Q):
                    #fstar[i,j,q] = (1 - Dt/thau)*f[i,j,q] + (Dt/thau)*feq[i,j,q]
                    fstar[i,j,q] = f[i,j,q] -  (f[i,j,q]-feq[i,j,q])/thau
                
    return fstar
def Streaming(f,Nx,Ny,Q,Bnd = None,with_boundaries = False):
    #Input 
    # f[i,j,q] ::: 
    # Nx :::
    # Ny :::
    #Output 
    # fstream[i,j,q] :::
    #Init variables 
    #Init variables
    if not with_boundaries :
        Bnd = np.zeros((Nx,Ny))
    fstream = np.zeros((Nx,Ny,Q))
    #print('Streaming step')
    for i in range(Nx): 
        for j in range(Ny):
            ip = (i + 1 + Nx) % Nx
            im = (i - 1 + Nx) % Nx
            jp = (j + 1 + Ny) % Ny
            jm = (j - 1 + Ny) % Ny
            
            #fstream[i,j,0] = f[i,j,0]
            #fstream[ip,j,1] = f[i,j,1]
            #fstream[i,jp,2] = f[i,j,2]
            #fstream[im,j,3] = f[i,j,3]
            #fstream[i,jm,4] = f[i,j,4]
            #fstream[ip,jp,5] = f[i,j,5]
            #fstream[im,jp,6] = f[i,j,6]
            #fstream[im,jm,7] = f[i,j,7]
            #fstream[ip,jm,8] = f[i,j,8]
            
            if Bnd[i,j] == 0 : 
                fstream[i,j,0] = f[i,j,0]
                fstream[i,j,1] = f[im,j,1]
                fstream[i,j,2] = f[im,jm,2]
                fstream[i,j,3] = f[i,jm,3]
                fstream[i,j,4] = f[ip,jm,4]
                fstream[i,j,5] = f[ip,j,5]
                fstream[i,j,6] = f[ip,jp,6]
                fstream[i,j,7] = f[i,jp,7]
                fstream[i,j,8] = f[im,jp,8]
    return fstream

def Bnd_HWBB(fold,f,label,normalx,normaly,Nx,Ny):
    fbnd = f 
    for i in range(Nx):
        for j in range(Ny):
            ip = (i + 1 + Nx) % Nx
            im = (i - 1 + Nx) % Nx
            jp = (j + 1 + Ny) % Ny
            jm = (j - 1 + Ny) % Ny
            if label[i,j] == 1 :
                #Bottom wall
                if normalx[i,j] == 0 and normaly[i,j] == 1 : 
                    fbnd[i,j,0] = fold[i,j,0]
                    fbnd[i,j,1] = fold[im,j,1]
                    fbnd[i,j,2] = fold[i,j,6] # HWBB
                    fbnd[i,j,3] = fold[i,j,7] # HWBB
                    fbnd[i,j,4] = fold[i,j,8] # HWBB
                    fbnd[i,j,5] = fold[ip,j,5]
                    fbnd[i,j,6] = fold[ip,jp,6]
                    fbnd[i,j,7] = fold[i,jp,7]
                    fbnd[i,j,8] = fold[im,jp,8]
                if normalx[i,j] == 0 and normaly[i,j] == -1 :
                    #Top Wall 
                    fbnd[i,j,0] = fold[i,j,0]
                    fbnd[i,j,1] = fold[im,j,1]
                    fbnd[i,j,2] = fold[im,jm,2]
                    fbnd[i,j,3] = fold[i,jm,3]
                    fbnd[i,j,4] = fold[ip,jm,4]
                    fbnd[i,j,5] = fold[ip,j,5]
                    fbnd[i,j,6] = fold[i,j,2] #
                    fbnd[i,j,7] = fold[i,j,3] #
                    fbnd[i,j,8] = fold[i,j,4] #   
                if normalx[i,j] == 1 and normaly[i,j] == 0 : 
                    #Left wall 
                    print('To be coded')
                    fbnd[i,j,0] = fold[i,j,0]
                    fbnd[i,j,1] = fold[im,j,1]
                    fbnd[i,j,2] = fold[im,jm,2]
                    fbnd[i,j,3] = fold[i,jm,3]
                    fbnd[i,j,4] = fold[ip,jm,4]
                    fbnd[i,j,5] = fold[ip,j,5]
                    fbnd[i,j,6] = fold[ip,jp,6]
                    fbnd[i,j,7] = fold[i,jp,7]
                    fbnd[i,j,8] = fold[im,jp,8]
                if normalx[i,j] == -1 and normaly[i,j] == 0 : 
                    #Right wall 
                    print('To be coded')
                    fbnd[i,j,0] = fold[i,j,0]
                    fbnd[i,j,1] = fold[im,j,1]
                    fbnd[i,j,2] = fold[im,jm,2]
                    fbnd[i,j,3] = fold[i,jm,3]
                    fbnd[i,j,4] = fold[ip,jm,4]
                    fbnd[i,j,5] = fold[ip,j,5]
                    fbnd[i,j,6] = fold[ip,jp,6]
                    fbnd[i,j,7] = fold[i,jp,7]
                    fbnd[i,j,8] = fold[im,jp,8]
    return fbnd

def WallBnd(fold,f,label,normalx,normaly,Nx,Ny): 
    fbnd = f
    for i in range(Nx):
        ip = (i + 1 + Nx) % Nx
        im = (i - 1 + Nx) % Nx
        jp = 1
        # bottom Wall 
        fbnd[i,0,0] = fold[0,0,0]
        fbnd[i,0,1] = fold[im,0,1]
        fbnd[i,0,2] = fold[i,0,6] #
        fbnd[i,0,3] = fold[i,0,7] #
        fbnd[i,0,4] = fold[i,0,8] #
        fbnd[i,0,5] = fold[ip,0,5]
        fbnd[i,0,6] = fold[ip,jp,6]
        fbnd[i,0,7] = fold[i,jp,7]
        fbnd[i,0,8] = fold[im,jp,8]
    for i in range(Nx):
        ip = (i + 1 + Nx) % Nx
        im = (i - 1 + Nx) % Nx
        jm = Ny-2
        # Top Wall 
        fbnd[i,j,0] = fold[i,j,0]
        fbnd[i,j,1] = fold[im,j,1]
        fbnd[i,j,2] = fold[im,jm,2]
        fbnd[i,j,3] = fold[i,jm,3]
        fbnd[i,j,4] = fold[ip,jm,4]
        fbnd[i,j,5] = fold[ip,j,5]
        fbnd[i,Ny-1,6] = fold[i,Ny-1,2]
        fbnd[i,Ny-1,7] = fold[i,Ny-1,3]
        fbnd[i,Ny-1,8] = fold[i,Ny-1,4]

    return fbnd


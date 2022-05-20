import numpy as np 

class Solver():
    def __init__(self):
        self.D = 0
        self.Q = 0
        self.ex = np.array([]) 
        self.ey = np.array([]) 
        self.ez = np.array([]) 
        self.w = np.array([]) 
        self.cs = 0. #(1/np.sqrt(3))

    def CalcFeq(self,rho,ux,uy,Nx,Ny,thau,srct =[0,0]):
        ## Input 
        # rho[i,j] :::
        # ux[i,j] :::
        # uy[i,j] :::
        # self.ex[q] :::
        # self.ey[q] :::
        # self.w[q] :::
        # self.cs ::: speed of sound 
        # Nx :::
        # Ny :::
        # self.Q :::
        ## Output 
        # feq[i,j,q] :::
        #Init variables
        
        ex = self.ex
        ey = self.ey
        ez = self.ez
        w = self.w
        Q = self.Q
        cs = self.cs

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

    def Collision(self,f,feq,thau,Nx,Ny,Bnd = None,with_boundaries = False):
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
        Q = self.Q
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

    def CalcMacro(self,f,Nx,Ny):
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
        ex = self.ex
        ey = self.ey
        Q = self.Q

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

class D2Q9(Solver):
    def __init__(self):
        self.D = 2
        self.Q = 9
        self.ex = np.array([0,1,1,0,-1,-1,-1,0,1])
        self.ey = np.array([0,0,1,1,1,0,-1,-1,-1])
        self.ez = np.array([]) 
        self.w =np.array([4/9,1/9,1/36,1/9,1/36,1/9,1/36,1/9,1/36])
        self.cs = (1/np.sqrt(3)) 

    def Streaming(self,f,Nx,Ny,Bnd = None,with_boundaries = False):
        #Input 
        # f[i,j,q] ::: 
        # Nx :::
        # Ny :::
        #Output 
        # fstream[i,j,q] :::
        #Init variables 
        Q = self.Q
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
    
    def WallBnd3(self,fold,f,label,normalx,normaly,Nx,Ny): 
        ex = self.ex
        ey = self.ey
        w = self.w
        Q = self.Q
        cs = self.cs
        
        fbnd = f
        for i in range(Nx):
            j = 0 
            ip = (i + 1 + Nx) % Nx
            im = (i - 1 + Nx) % Nx
            jp = 1
            # bottom Wall 
            fbnd[i,j,0] = fold[i,j,0]
            fbnd[i,j,1] = fold[im,j,1]
            fbnd[i,j,2] = fold[i,j,6] # HWBB
            fbnd[i,j,3] = fold[i,j,7] # HWBB
            fbnd[i,j,4] = fold[i,j,8] # HWBB
            fbnd[i,j,5] = fold[ip,j,5]
            fbnd[i,j,6] = fold[ip,jp,6]
            fbnd[i,j,7] = fold[i,jp,7]
            fbnd[i,j,8] = fold[im,jp,8]
        for i in range(Nx):
            j = Ny - 1 
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
        
        j = 1 
        while j < Ny-1:
        #for j in range(Ny):
            i = 0 
            ip = (i + 1 + Nx) % Nx # 1
            #im = (i - 1 + Nx) % Nx # Doesn't exists 
            jp = (j + 1 + Ny) % Ny
            jm = (j - 1 + Ny) % Ny
            
            # Left Wall 
            fbnd[i,j,0] = fold[i,j,0] 
            fbnd[i,j,3] = fold[i,jm,3]
            fbnd[i,j,4] = fold[ip,jm,4]
            fbnd[i,j,5] = fold[ip,j,5]
            fbnd[i,j,6] = fold[ip,jp,6]
            fbnd[i,j,7] = fold[i,jp,7]
            
            cs_squared = cs * cs 
            Uwx = 0.1 * cs 
            Uwy = 0.
            dens = 1.
            fbnd[i,j,1] = fold[i,j,5] -  w[5]*dens*(ex[5]*Uwx + ey[5]*Uwy)/cs_squared#+ Uw
            fbnd[i,j,2] = fold[i,j,6] -  w[6]*dens*(ex[6]*Uwx + ey[6]*Uwy)/cs_squared#-  (1/(6*c))*Uw#
            fbnd[i,j,8] = fold[i,j,4] -  w[4]*dens*(ex[4]*Uwx + ey[4]*Uwy)/cs_squared#+  (1/(6*c))*Uw#
            j = j + 1
        
        j = 1
        while j < Ny - 1 : 
            i = Nx-1
            #ip = (i + 1 + Nx) % Nx # Doesn't exists
            im = (i - 1 + Nx) % Nx # Nx-2
            jp = (j + 1 + Ny) % Ny
            jm = (j - 1 + Ny) % Ny
            #right wall
            for q in range(Q):
                fbnd[i,j,q] = fold[im,j,q]
            j = j + 1
        return fbnd

    def Bnd_HWBB(self,fold,f,label,normalx,normaly,Nx,Ny,Inletx,Inlety):
        fbnd = f 
        w = self.w
        ex = self.ex
        ey = self.ey
        cs = self.cs
        Q = self.Q
        cs_squared = cs * cs 
        for i in range(Nx):
            for j in range(Ny):
                ip = (i + 1 + Nx) % Nx
                im = (i - 1 + Nx) % Nx
                jp = (j + 1 + Ny) % Ny
                jm = (j - 1 + Ny) % Ny
                if label[i,j] == 1 :
                    #Bottom wall
                    if normalx[i,j] == 0 and normaly[i,j] == 1 : 
                        Uwx = Inletx[i,j] #0.1*cs #Inletx[i,j]
                        Uwy = Inlety[i,j] #0. #Inlety[i,j]
                        dens = 1.
                        fbnd[i,j,0] = fold[i,j,0]
                        fbnd[i,j,1] = fold[im,j,1]
                        fbnd[i,j,2] = fold[i,j,6] - w[6]*dens*(ex[6]*Uwx + ey[6]*Uwy)/cs_squared# HWBB
                        fbnd[i,j,3] = fold[i,j,7] - w[7]*dens*(ex[7]*Uwx + ey[7]*Uwy)/cs_squared# HWBB
                        fbnd[i,j,4] = fold[i,j,8] - w[8]*dens*(ex[8]*Uwx + ey[8]*Uwy)/cs_squared# HWBB
                        fbnd[i,j,5] = fold[ip,j,5]
                        fbnd[i,j,6] = fold[ip,jp,6]
                        fbnd[i,j,7] = fold[i,jp,7]
                        fbnd[i,j,8] = fold[im,jp,8]

                    if normalx[i,j] == 0 and normaly[i,j] == -1 :
                        Uwx = Inletx[i,j] #0.1*cs #Inletx[i,j]
                        Uwy = Inlety[i,j] #0. #Inlety[i,j]
                        dens = 1.
                        #Top Wall 
                        fbnd[i,j,0] = fold[i,j,0]
                        fbnd[i,j,1] = fold[im,j,1]
                        fbnd[i,j,2] = fold[im,jm,2]
                        fbnd[i,j,3] = fold[i,jm,3]
                        fbnd[i,j,4] = fold[ip,jm,4]
                        fbnd[i,j,5] = fold[ip,j,5]
                        fbnd[i,j,6] = fold[i,j,2] - w[2]*dens*(ex[2]*Uwx + ey[2]*Uwy)/cs_squared# HWBB # HWBB
                        fbnd[i,j,7] = fold[i,j,3] - w[3]*dens*(ex[3]*Uwx + ey[3]*Uwy)/cs_squared# HWBB# HWBB
                        fbnd[i,j,8] = fold[i,j,4] - w[4]*dens*(ex[4]*Uwx + ey[4]*Uwy)/cs_squared# HWBB# HWBB

                    if normalx[i,j] == 1 and normaly[i,j] == 0 : 
                        #Left wall 
                        #print('To be coded')
                        #fbnd[i,j,0] = fold[i,j,0]
                        #fbnd[i,j,1] = fold[im,j,1]
                        #fbnd[i,j,2] = fold[im,jm,2]
                        #fbnd[i,j,3] = fold[i,jm,3]
                        #fbnd[i,j,4] = fold[ip,jm,4]
                        #fbnd[i,j,5] = fold[ip,j,5]
                        #fbnd[i,j,6] = fold[ip,jp,6]
                        #fbnd[i,j,7] = fold[i,jp,7]
                        #fbnd[i,j,8] = fold[im,jp,8]

                        Uwx = Inletx[i,j] #0.1*cs #Inletx[i,j]
                        Uwy = Inlety[i,j] #0. #Inlety[i,j]
                        dens = 1.

                        fbnd[i,j,0] = fold[i,j,0] 
                        fbnd[i,j,3] = fold[i,jm,3]
                        fbnd[i,j,4] = fold[ip,jm,4]
                        fbnd[i,j,5] = fold[ip,j,5]
                        fbnd[i,j,6] = fold[ip,jp,6]
                        fbnd[i,j,7] = fold[i,jp,7]
                        fbnd[i,j,1] = fold[i,j,5] -  w[5]*dens*(ex[5]*Uwx + ey[5]*Uwy)/cs_squared# HWBB
                        fbnd[i,j,2] = fold[i,j,6] -  w[6]*dens*(ex[6]*Uwx + ey[6]*Uwy)/cs_squared# HWBB
                        fbnd[i,j,8] = fold[i,j,4] -  w[4]*dens*(ex[4]*Uwx + ey[4]*Uwy)/cs_squared# HWBB
                        
                    if normalx[i,j] == -1 and normaly[i,j] == 0 : 
                        #Right wall 
                        #print('To be coded')
                        #fbnd[i,j,0] = fold[i,j,0]
                        #fbnd[i,j,1] = fold[im,j,1]
                        #fbnd[i,j,2] = fold[im,jm,2]
                        #fbnd[i,j,3] = fold[i,jm,3]
                        #fbnd[i,j,4] = fold[ip,jm,4]
                        #fbnd[i,j,5] = fold[ip,j,5]
                        #fbnd[i,j,6] = fold[ip,jp,6]
                        #fbnd[i,j,7] = fold[i,jp,7]
                        #fbnd[i,j,8] = fold[im,jp,8]

                        Uwx = Inletx[i,j] #0.1*cs #Inletx[i,j]
                        Uwy = Inlety[i,j] #0. #Inlety[i,j]
                        dens = 1.

                        fbnd[i,j,0] = fold[i,j,0] 
                        fbnd[i,j,1] = fold[im,j,5] 
                        fbnd[i,j,2] = fold[im,jm,6] 
                        fbnd[i,j,3] = fold[i,jm,3] 
                        fbnd[i,j,7] = fold[i,jp,7]
                        fbnd[i,j,8] = fold[im,jp,4] 
                        fbnd[i,j,4] = fold[i,j,8] -  w[8]*dens*(ex[8]*Uwx + ey[8]*Uwy)/cs_squared#
                        fbnd[i,j,5] = fold[i,j,1] -  w[1]*dens*(ex[1]*Uwx + ey[1]*Uwy)/cs_squared#
                        fbnd[i,j,6] = fold[i,j,2] -  w[2]*dens*(ex[2]*Uwx + ey[2]*Uwy)/cs_squared#
        j = 1
        while j < Ny - 1 : 
            i = Nx-1
            if label[i,j] == 2:
                #ip = (i + 1 + Nx) % Nx # Doesn't exists
                im = (i - 1 + Nx) % Nx # Nx-2
                jp = (j + 1 + Ny) % Ny
                jm = (j - 1 + Ny) % Ny
                #right wall
                for q in range(Q):
                    fbnd[i,j,q] = fold[im,j,q]
            j = j + 1

        return fbnd

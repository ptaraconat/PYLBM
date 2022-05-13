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
        feq = np.zeros((Nx,Ny,Q))
        ex = self.ex
        ey = self.ey
        ez = self.ez
        w = self.w
        Q = self.Q
        cs = self.cs

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

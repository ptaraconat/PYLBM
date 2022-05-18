import matplotlib.pyplot as plt 
import sys as sys 
sys.path.append('../src/')
from LaticeNumerics import * 

##Physical parameters 
nu = 1
Pref = 101325
Tref = 300
Rgas = 287.15

##Numerical parameters 
Solver = D2Q9() 
Dt = 0.000001
Dx = 0.015
# Calc tau
Nulb = (nu*Dt)/(Dx * Dx)
thau = (Nulb/(Solver.cs*Solver.cs)) + 0.5

##Create Fluid domain and Grid 
Lx = 0.9
Ly = 0.3
Nx = np.int((Lx/Dx) + 1)
Ny = np.int((Ly/Dx) + 1)
x = np.linspace(-Lx/2,Lx/2,Nx)
y = np.linspace(-Ly/2,Ly/2,Ny)
xv, yv = np.meshgrid(x, y,indexing = 'ij') 

##Boundary conditions 
label = np.zeros((Nx,Ny))
label[:,0] = 1 
label[:,Ny-1] = 1 
label[0,:] = 1 
label[Nx-1,:] = 1 
#
normaly = np.zeros((Nx,Ny))
normalx = np.zeros((Nx,Ny))
normaly[:,0] = 1
normalx[:,0] = 0
normaly[:,Ny-1] = -1
normalx[:,Ny-1] = 0

##Initial condition 
ux = np.zeros((Nx,Ny))
uy = np.zeros((Nx,Ny))
P = np.ones((Nx,Ny)) * Pref
rho = P/(Rgas*Tref) 

# Scale Macro var 
rho = rho *((Rgas*Tref)/Pref)

#Temporal loop
f = Solver.CalcFeq(rho,ux,uy,Nx,Ny,thau,srct =[0,0])
for i in range(1001):
    feq = Solver.CalcFeq(rho,ux,uy,Nx,Ny,thau,srct =[0,0])
    fstar = Solver.Collision(f,feq,thau,Nx,Ny,Bnd = label,with_boundaries = True)
    fold = fstar 
    fstream = Solver.Streaming(fstar,Nx,Ny,Bnd = label,with_boundaries = True)
    f = fstream
    #fbnd = WallBnd(fold,f,label,normalx,normaly,Nx,Ny)
    fbnd = Solver.WallBnd3(fold,f,label,normalx,normaly,Nx,Ny)
    f = fbnd
    ux,uy,rho = Solver.CalcMacro(fstream,Nx,Ny)
    
    if i%10 == 0 :
        print(i)
    if i%100 == 0:
        plt.contourf(xv,yv,ux,cmap = 'jet')
        plt.colorbar()
        plt.savefig('ChannelV2'+str(i)+'png',format = 'png')
        plt.show()
        plt.close()
        


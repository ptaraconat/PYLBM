import matplotlib.pyplot as plt 
import sys as sys 
sys.path.append('../src/')
from LaticeNumerics import * 

##Physical parameters 
nu = 0.01
Pref = 101325
Tref = 300
Rgas = 287.15
vforce = 0

##Numerical parameters 
Solver = D2Q9() 
Dt = 0.00001
Dx = 0.001
# Calc tau
Nulb = (nu*Dt)/(Dx * Dx)
thau = (Nulb/(Solver.cs*Solver.cs)) + 0.5

##Create Fluid domain and Grid 
Lx = 0.1
Ly = 0.1
Nx = np.int((Lx/Dx) + 1)
Ny = np.int((Ly/Dx) + 1)
x = np.linspace(-Lx/2,Lx/2,Nx)
y = np.linspace(-Ly/2,Ly/2,Ny)
xv, yv = np.meshgrid(x, y,indexing = 'ij') 

##Boundary conditions 
label = np.zeros((Nx,Ny))
normaly = np.zeros((Nx,Ny))
normalx = np.zeros((Nx,Ny))
Uwx = np.zeros((Nx,Ny))
Uwy = np.zeros((Nx,Ny))

##Initial condition 
#Initial condition 
ux = np.zeros((Nx,Ny))
uy = np.zeros((Nx,Ny))
rho = np.zeros((Nx,Ny))
P = np.zeros((Nx,Ny))
for i in range(Nx):
    for j in range(Ny):
        R = 0.004/Dx
        xCentered =  i-(Nx-1.)/2.#xv[i,j] 
        yCentered =  j-(Ny-1.)/2.#yv[i,j]
        P[i,j] = Pref+ 0.1 * Pref * np.exp(-(xCentered*xCentered+yCentered*yCentered)/(4*R*R))
rho = P/(Rgas*Tref) 

# Scale Macro var 
rho = rho *((Rgas*Tref)/Pref)
vforce = vforce * ((Dt*Dt)/Dx)

#Temporal loop
f = Solver.CalcFeq(rho,ux,uy,Nx,Ny,thau,srct =[vforce,0])
for i in range(101):
    feq = Solver.CalcFeq(rho,ux,uy,Nx,Ny,thau,srct =[vforce,0])
    fstar = Solver.Collision(f,feq,thau,Nx,Ny,Bnd = label,with_boundaries = True)
    fold = fstar 
    fstream = Solver.Streaming(fstar,Nx,Ny,Bnd = label,with_boundaries = True)
    f = fstream
    #fbnd = Solver.WallBnd4(fold,f,label,normalx,normaly,Nx,Ny)
    fbnd = Solver.Bnd_HWBB(fold, f, label, normalx, normaly, Nx, Ny,Uwx,Uwy)
    f = fbnd
    ux,uy,rho = Solver.CalcMacro(fstream,Nx,Ny)
    
    if i%1 == 0 :
        print(i)
    # Post Treatement 
    if i%10 == 0:
        print('dumping')
        fig = plt.figure(figsize = (12,4))
        ax1 = fig.add_subplot(131)
        ax2 = fig.add_subplot(132)
        ax3 = fig.add_subplot(133)
        z1_plot=ax1.contourf(xv, yv, rho,cmap = 'jet')
        ax2.contourf(xv, yv, ux,cmap = 'jet')
        ax3.contourf(xv, yv, uy,cmap = 'jet')
        #fig.colorbar(z1_plot,cax=ax1)
        #plt.show()
        plt.savefig('lbmpulse'+str(i)+'png',format = 'png')
        plt.close()
        


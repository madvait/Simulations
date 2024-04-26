# here, you may find my not so efficient simulation of lid driven cavity flow which I implemented for my understanding of streaming and moving wall boundary conditions


import numpy as np
import matplotlib.pyplot as plt
from lbm import *
from tqdm import tqdm

#problem constant
nx,ny = 41,40

lid_velocity = 0.4
iterationtime = 10000
Re = 500
# lattice constants
nu = (lid_velocity*ny)/Re
tau = 3*nu+0.5
#tau = 0.525
w = np.array([4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36])
c = np.array([[0,0],[1, 0],[0, 1],[-1, 0],[0, -1],[1, 1],[-1, 1],[-1, -1],[1, -1]])
reverse_indices = np.array([0,3,4,1,2,7,8,5,6])

# initialize some constants:
rho = np.ones([nx,ny])
u = np.zeros([nx,ny,2])
f = calc_feq_tensor(u,rho)

def update(f_prev):

    rho = calc_density(f_prev)
    rho[:,0] = 1
    rho[(0,-1),:] = 1
    u = calc_velocity_vectorspace(f_prev,rho)
    u[:,0,:] = 0
    u[(0,-1),:,:]=0
    f_eq = calc_feq_tensor(u,rho)
    f_prev[:,0,:] = f_eq[:,0,:]
    f_prev[(0,-1),:,:]= f_eq[(0,-1),:,:]

    f_star = f_prev - (1/tau)*(f_prev-f_eq)

    f_temporary = np.copy(f_star) #a temporary matrix to store streaming data
    #stream at the fluid range
    for i in range(1,nx-1):
        for j in range(1,ny):
            for a in range(9):
                ia = i+c[a,0]
                ja = j+c[a,1]
                if ja>ny-1: ja =0 #garbage value, to be collected
                f_temporary[ia,ja,a] = f_star[i,j,a] #this should stream all the f not a part of the wall

                #note: the values in the corner are not never streamed and should be replaced with the equilibrium value for proper calculations of velocity and density
                f_temporary[1,1,8] = f_eq[1,1,8]
                f_temporary[1,1,6] = f_eq[1,1,6]

                f_temporary[-2,1,5] = f_eq[-2,1,5]
                f_temporary[-2,1,7] = f_eq[-2,1,7]

                f_temporary[1,-1,5] = f_eq[1,-1,5]
                f_temporary[1,-1,7] = f_eq[1,-1,7]

                f_temporary[-2,-1,8] = f_eq[-2,-1,8]
                f_temporary[-2,-1,6] = f_eq[-2,-1,6]
                
                
        #for the boundaries of i=0,i=nx-1 and j = 0:
        #manually apply the boundary conditions for i->-1,1 and j =1
    f_temporary[1:-1,1,5] = f_star[1:-1,1,7]
    f_temporary[1:-1,1,2] = f_star[1:-1,1,4]
    f_temporary[1:-1,1,6] = f_star[1:-1,1,8]

    # for j = 1,: at the left side i = 1
    f_temporary[1,1:,1] = f_star[1,1:,3]
    f_temporary[1,1:,5] = f_star[1,1:,7]
    f_temporary[1,1:,8] = f_star[1,1:,6]

    #for j = 1,: at the right side i =-2
    f_temporary[-2,1:,3] = f_star[-2,1:,1]
    f_temporary[-2,1:,7] = f_star[-2,1:,5]
    f_temporary[-2,1:,6] = f_star[-2,1:,8]
    
    f_stream = np.copy(f_temporary)
    # calculate the boundary density for appying the moving wall boundary condition
    rho[:,-1] = f_stream[:,-1,0]+f_stream[:,-1,1]+f_stream[:,-1,3]+2*(f_stream[:,-1,2]+f_stream[:,-1,6]+f_stream[:,-1,5])
    # do check if it is necessary to get it at the boundary i = 1 to -1

    f_stream[2:-2,-1,4] = f_stream[2:-2,-1,2]
    f_stream[2:-2,-1,7] = f_stream[2:-2,-1,5]- ((rho[2:-2,-1]*lid_velocity)/6)
    f_stream[2:-2,-1,8] = f_stream[2:-2,-1,6]+ ((rho[2:-2,-1]*lid_velocity)/6)

    return f_stream


x = np.arange(nx)
y = np.arange(ny)
x,y = np.meshgrid(x, y)

plt.figure(figsize=(5,4), dpi = 100)

for time in tqdm(range(iterationtime)):
    fnext = update(f)
    f = np.copy(fnext)
    rho = calc_density(fnext)
    u = calc_velocity_vectorspace(fnext,rho)
    if time%3==0 and time>200:
        velocity_magnitude = np.linalg.norm(u, axis=-1)
        plt.subplot(111)
        plt.contourf(x,y,np.transpose(velocity_magnitude), cmap='magma')
        plt.colorbar(label='Velocity Magnitude')
        plt.title('Velocity Magnitude Map')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.draw()
        plt.pause(0.0001)
        plt.clf()
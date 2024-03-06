## Cellular Automata Simulations of Conway's game of life
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

Nx= 120
Ny= 120

#defining a function that outputs the next generation of the game

def nextgen(gen1,Nx,Ny):
    gen2 = np.copy(gen1)
    living = np.zeros((Nx,Ny))
    for i in range(Nx):
        for j in range(Ny):
            if i!=Nx-1 and j!= Ny-1:
                #calculate the number of living neighbours
                living[i,j] = gen1[i-1,j-1]+gen1[i-1,j]+gen1[i-1,j+1]+gen1[i,j-1]+gen1[i,j+1]+gen1[i+1,j-1]+gen1[i+1,j]+gen1[i+1,j+1]
            elif i==Nx-1 and j!=Ny-1:
                living[i,j] = gen1[i-1,j-1]+gen1[i-1,j]+gen1[i-1,j+1]+gen1[i,j-1]+gen1[i,j+1]+gen1[0,j-1]+gen1[0,j]+gen1[0,j+1]
            elif i!=Nx-1 and j==Ny-1:
                living[i,j] = gen1[i-1,j-1]+gen1[i-1,j]+gen1[i-1,0]+gen1[i,j-1]+gen1[i,0]+gen1[i+1,j-1]+gen1[i+1,j]+gen1[i+1,0]
            elif i==Nx-1 and j==Ny-1:
                living[i,j] = gen1[i-1,j-1]+gen1[i-1,j]+gen1[i-1,0]+gen1[i,j-1]+gen1[i,0]+gen1[0,j-1]+gen1[0,j]+gen1[0,0]

                #rules for propagation
    for i in range(Nx):
        for j in range(Ny):
            if living[i,j]<2 and gen1[i,j]==1:
                gen2[i,j]=0
            elif (living[i,j]==2 or living[i,j]==3) and gen1[i,j]==1:
                gen2[i,j]=1
            elif living[i,j]>3 and gen1[i,j]==1:
                gen2[i,j]=0
            elif living[i,j]==3 and gen1[i,j]==0:
                gen2[i,j]=1

    return gen2




# Initialize your variables
gen1 = np.random.randint(2, size=(Nx, Ny))


# Function to update the plot for each frame
def update(frame):
    global gen1

    gen1 = nextgen(gen1,np.shape(gen1)[0],np.shape(gen1)[1])

    im.set_array(gen1)  # Update the image data

    return im,

# Create a figure and axis
fig, ax = plt.subplots()
im = ax.imshow(gen1, cmap='gray')

# Create a FuncAnimation object to generate the animation
ani = FuncAnimation(fig, update, frames= True, interval=20)  # Change frames and interval as needed
plt.axis('off')
plt.show()


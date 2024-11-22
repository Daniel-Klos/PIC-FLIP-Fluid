# PIC-FLIP-Fluid
A multithreaded particle in cell (PIC) fluid simulation combined with a fluid implicit particle (FLIP) simulation. 

There are two different types of collision detection available, one which uses constant memory for every grid size, and another which increases the amount of memory used when the grid size is increased. For the most speed, use non-constant memory collision detection, as it is multithreaded, but if changing the grid size (not the amount of particles) during runtime with no change in memory is a must then go for constant memory.
In the simulate() method, just comment out the method of collision detection that you don't want to use.

Controls:
-   Press 1 for an interactive rigid object and 2 for a force object.
-   Press s to bias the simulation towards pic, and b to bias the simulation towards flip.
-   Press r to decrease the force object radius, and t to increase
-   Press w to decrease the vorticity, and e to increase it
-   Press d to decrease the divergence modifier, c to increase it
-   Press f to switch between velocity rendering and diffusion rendering
-   Press q to exit the simulation

This implementation includes vorticity confinement.

The simulation is able to run 30k particles at around 80 fps, as you can see in the pictures below.

![Screenshot 2024-11-11 000055](https://github.com/user-attachments/assets/b5d76477-b647-44a8-98c6-88b33f32bd13)

![Screenshot 2024-11-11 000239](https://github.com/user-attachments/assets/7cdf9d67-af4b-419b-a21e-06990b2f07d0)

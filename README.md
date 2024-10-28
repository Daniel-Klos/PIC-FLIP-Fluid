# PIC-FLIP-Fluid
A particle in cell (PIC) fluid simulation combined with a fluid implicit particle (FLIP) simulation. 
I am currently working on multithreading collision detection which will hopefully speed up this project.

There are two different types of collision detection available, one which uses constant memory for every grid size, and another which increases the amount of memory used when the grid size is increased. For the most speed, use non-constant memory collision detection, but if changing the grid size (not the amount of particles) is a must then go for constant memory.
In the simulate() method, just comment out the method of collision detection that you don't want to use. Do not mismatch, you must initialize and query a constant memory spatial hash if that is what you are going to use. 

The simulation is able to run 10k particles at a smooth ~85 fps, as you can see in the picture below.

![Screenshot 2024-10-27 200111](https://github.com/user-attachments/assets/c5c37cd4-b8ee-4933-a0a7-10611eb0007e)

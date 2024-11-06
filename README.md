# PIC-FLIP-Fluid
A particle in cell (PIC) fluid simulation combined with a fluid implicit particle (FLIP) simulation. 
I am currently working on improving the multithreaded collision detection which will hopefully speed up this project.

There are two different types of collision detection available, one which uses constant memory for every grid size, and another which increases the amount of memory used when the grid size is increased. For the most speed, use non-constant memory collision detection, but if changing the grid size (not the amount of particles) during runtime with no change in memory is a must then go for constant memory.
In the simulate() method, just comment out the method of collision detection that you don't want to use. Do not mismatch, you must initialize and query a constant memory spatial hash if that is what you are going to use. 

The simulation is able to run 30k particles at close to 60 fps, as you can see in the picture below.

![Screenshot 2024-10-30 010439](https://github.com/user-attachments/assets/e057f33d-1b80-4e07-a6d5-76bfe07f58f9)
![Screenshot 2024-11-06 011221](https://github.com/user-attachments/assets/24636d32-62e4-401f-91ae-5ae7759be2c0)

# PIC-FLIP-Fluid
A 99% multithreaded particle in cell (PIC) fluid simulation combined with a fluid implicit particle (FLIP) simulation.

Controls:
-   Press 1 for an interactive rigid object, 2 for a force object, 3 for an object that generates/removes particles, and 4 to draw/erase solids (left & right click). Use the mouse to change the size of objects.
-   Press s to bias the simulation towards pic, and b to bias the simulation towards flip.
-   Use the mouse wheel to increase/decrease the radius of objects
-   Press w to decrease the vorticity confinement, and e to increase it
-   Press d to decrease the divergence modifier, c to increase it
-   Press a to switch between diffusion rendering, velocity rendering, vorticity rendering, temperature rendering, and divergence rendering
    - While in temperature rendering mode, press f to turn on fire mode
-   Press q to exit the simulation

The simulation is able to run around 200k particles at 60 fps

![Screenshot 2025-05-31 144513](https://github.com/user-attachments/assets/9998471c-b16c-45f5-b34b-5f588ca8f96a)

![Screenshot 2025-02-23 142512](https://github.com/user-attachments/assets/7eb92834-c2e5-4f8c-90be-1b4fce782517)

![Screenshot 2025-02-23 142916](https://github.com/user-attachments/assets/5981318f-5ac5-4cc7-aacd-515fefe743ca)

![Screenshot 2025-03-26 005615](https://github.com/user-attachments/assets/988be616-f1b5-483a-9f0e-76a55853a383)

Improvements to be made:
-  implement Multigrid, and then MGPCG
    - Efficient Red Black SOR/GS has been implemented as well as general Preconditioned Conjugate Gradient code, now need to implement efficient Dampened Jacobi
-  DDA ray casting & making sure that particle-particle collisions don't push particles into obstacles
-  Actual fluid rendering, in addition to the current particle view
-  Implicit Density Projection
-  Currently using no-slip boundary conditions. Implement slip boundary conditions from that Muller paper
-  Use thread buffers more (to grid most importantly)

# PIC-FLIP-Fluid
A multithreaded particle in cell (PIC) fluid simulation combined with a fluid implicit particle (FLIP) simulation. 

Controls:
-   Press 1 for an interactive rigid object, 2 for a force object, 3 for an object that generates/removes particles, and 4 to draw/erase solids (left & right click)
-   Press s to bias the simulation towards pic, and b to bias the simulation towards flip.
-   Use the mouse wheel to increase/decrease the radius of objects
-   Press w to decrease the vorticity confinement, and e to increase it
-   Press d to decrease the divergence modifier, c to increase it
-   Press a to switch between diffusion rendering, velocity rendering, vorticity rendering and temperature rendering
    - While in temperature rendering mode, press f to turn on fire mode
-   Press q to exit the simulation

This implementation also includes vorticity confinement.

The simulation is able to run 35k particles at around 80 fps, as you can see in the pictures below.

![Screenshot 2025-02-23 142512](https://github.com/user-attachments/assets/7eb92834-c2e5-4f8c-90be-1b4fce782517)

![Screenshot 2025-02-23 142722](https://github.com/user-attachments/assets/b5f6bdc6-fa71-4b44-9ea7-9f305d571fb6)

![Screenshot 2025-02-23 142916](https://github.com/user-attachments/assets/5981318f-5ac5-4cc7-aacd-515fefe743ca)

![Screenshot 2025-03-26 005615](https://github.com/user-attachments/assets/988be616-f1b5-483a-9f0e-76a55853a383)

Improvements to be made:
-  Finish preconditioned conjugate gradient code and optimize it
-  complete multithreading
-  better collision resolution with scene objects
-  more accurate volume correction

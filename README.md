# High-resolution schemes
- Contain both 5th-order-WENO (Weighted Essentially Non-Oscillatory) and MUSCL (Monotonic Upstream-centered Scheme for Conservation Laws) schemes
- Hydrodynamic or Magnetohydrodynamic (MHD) solver. 
- Second-order accurate Runge-Kutta time evolution
- Viscous effects up to 4th-order accuracy

1. 1st-order Godunov
2. 2nd-order Monotonized Central (MC) 
3. 5th-order WENO (JS)


Example simulation on shock-bubble interaction (400 by 200) using 5th-order WENO with symmetry boundary condition: <br />
https://www.youtube.com/watch?v=HPZ6TjSnV1I&t=1s&ab_channel=Spaceduck496 <br />
(a) For inviscid flow at t=0.15s,
![den_n_p](https://user-images.githubusercontent.com/64028216/216842603-724e618c-1cf5-49af-851f-1bd4fae4223e.png)
![vel](https://user-images.githubusercontent.com/64028216/216842635-2e14e490-845e-483b-a142-1c99c492d6bc.png)

(b) For viscous flow at t=1.0s, the temperature field for viscous shock tube is,
![Screenshot from 2023-02-20 15-34-48](https://user-images.githubusercontent.com/64028216/220245447-f0ca847b-8a9e-4c2a-8019-362514f03d8a.png)



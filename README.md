# High-resolution schemes
- Contain both 5th-order-WENO (Weighted Essentially Non-Oscillatory) and MUSCL (Monotonic Upstream-centered Scheme for Conservation Laws) schemes
- Hydrodynamic or Magnetohydrodynamic (MHD) solver. Make use of dynamic memory management to minimize the amount of memory required.

1. 1st-order Godunov
2. 2nd-order Monotonized Central (MC) 
3. 5th-order WENO (JS)


Example simulation on shock-bubble interaction (400 by 200) using 5th-order WENO: <br />
https://www.youtube.com/watch?v=HPZ6TjSnV1I&t=1s&ab_channel=Spaceduck496 <br />
At t=0.15s,
![den_n_p](https://user-images.githubusercontent.com/64028216/216842603-724e618c-1cf5-49af-851f-1bd4fae4223e.png)
![vel](https://user-images.githubusercontent.com/64028216/216842635-2e14e490-845e-483b-a142-1c99c492d6bc.png)


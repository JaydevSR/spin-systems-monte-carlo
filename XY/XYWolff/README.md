# Wolff Algorithm for the XY Model

The lattice is composed of floating point numbers from 0 to 1. Each numeber represents the angle of spin w.r.t the x-axis normalized by `2PI`.

1. Choose a seed at random from the lattice. Also choose a random angle (unit vector) in the plane to project the spins for rotation (say `r`).

2. Flip the seed spin using the transformation: `s --> 2*r - s + 0.5` (in terms of angle / 2PI).

3. Look for the nearest neighbors of the seed spin say `n`, and add them to the cluster with probability `P_add = 1 - exp(-2 * cos(r - s) * cos(r - n))`. 

4. Repeat step `(2)` and `(3)` for all the spins that are added to the cluster. This step is repeated till there are no more spins to consider. If some spin is already added then it is not considered again. If some spin is already considered but not added then it can be considered again as the neighbour of another spin.
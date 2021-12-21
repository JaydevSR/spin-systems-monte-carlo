# Wolff Algorithm for the Ising Model

1. Choose a seed at random from the lattice.

2. Look at the neighbors of the spin, if any neighboring spin is pointing in the same direction the add it to the cluster with probability `P_add = 1-e^(-2JT)`.

3. Repeat step `(2)` for all the spins that are added to the cluster. This step is repeated till there are no more spins to consider. If some spin is already added then it is not considered again. If some spin is already considered but not added then it can be considered again as the neighbour of another spin.

4. Flip the cluster.
# Spin-System-Monte-Carlo

Code for preparing spin lattices (Ising and XY) using Monte-Carlo-Markov-Chains (Metropolis and Wolff algoritm) and making measurements on them.

## Metropolis Algorithm for Ising Model

1. We take a lattice with ***periodic boundary conditions***  and `B = 0`. Usually the initial state is choosen to be one at temperature `T = 0` where all the spins are alligned in same direction or one at temperature `T = \infty` where they are aligned randomly.

2. Choose a temperature to perform the simulation for. If the simulation needs to be performed for more then one temperature values then start at some value and use its final equilibrium state as the initial state of next temperature value. We take `\beta = 1 / T` i.e. `k = 1`.

3. The difference between the energies of two consecutive states such that a spin at site `k` is flipped to get the next state is given by (can be easily calculated by substitution and cancellation):
    ```
        delta_E = E_\nu - E_\mu = 2 * J * s_k^\mu \sum_{i n.n. to k} s_i^\mu
    ```
    
4. Now if `delta_E <= 0` we flip the spin otherwise we flip the spin with probability `e^{-delta_E / T}`. For this we pick a number `r` such that `0 <= r < 1` and flip the spin if `r < e^{- delta_E / T}`.

5. Then we run the simulation for several sweep (`1 sweep = N flip attempts`) and once the system reaches equilibrium we calculate the desired quantity. In order to determine whether the system has reached equilibrium or not we can look at the change in magnetization or internal energy with each sweep as both the quantities become stable on equilibrium.


## Wolff Algorithm for the Ising Model

1. Choose a seed at random from the lattice.

2. Look at the neighbors of the spin, if any neighboring spin is pointing in the same direction the add it to the cluster with probability `P_add = 1-e^(-2JT)`.

3. Repeat step `(2)` for all the spins that are added to the cluster. This step is repeated till there are no more spins to consider. If some spin is already added then it is not considered again. If some spin is already considered but not added then it can be considered again as the neighbour of another spin.

4. Flip the cluster.

## Wolff Algorithm for the XY Model

The lattice is composed of floating point numbers from 0 to 1. Each numeber represents the angle of spin w.r.t the x-axis normalized by `2PI`.

1. Choose a seed at random from the lattice. Also choose a random angle (unit vector) in the plane to project the spins for rotation (say `r`).

2. Flip the seed spin using the transformation: `s --> 2*r - s + 0.5` (in terms of angle / 2PI).

3. Look for the nearest neighbors of the seed spin say `n`, and add them to the cluster with probability `P_add = 1 - exp(-2 * cos(r - s) * cos(r - n))`. 

4. Repeat step `(2)` and `(3)` for all the spins that are added to the cluster. This step is repeated till there are no more spins to consider. If some spin is already added then it is not considered again. If some spin is already considered but not added then it can be considered again as the neighbour of another spin.
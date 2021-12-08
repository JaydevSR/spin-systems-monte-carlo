# Metropolis Algorithm for Ising Model

- We take a lattice with ***periodic boundary conditions***  and $B=0$. Usually the initial state is choosen to be one at temperature $T=0$ where all the spins are alligned in same direction or one at temperature $T=\infty$ where they are aligned randomly.
- Choose a temperature to perform the simulation for. If the simulation needs to be performed for more then one temperature values then start at some value and use its final equilibrium state as the initial state of next temperature value. We take $\beta = 1 / T$ i.e. $k=1$.
- The difference between the energies of two consecutive states such that a spin at site $k$ is flipped to get the next state is given by (can be easily calculated by substitution and cancellation):
    
    $$E_\nu - E_\mu = 2Js_k^\mu \sum_{i \ \text{n.n. to } k} s_i^\mu$$
    
- Now if $E_\nu - E_\mu \le$ 0 we flip the spin other wise we flip the spin with probability $e^{-(E_\nu-E_\mu) / T}$. For this we pick a number $r$ such that $0 \le r \lt 1$ and flip the spin if $r \lt e^{-(E_\nu-E_\mu) / T}$.

- Then we run the simulation for several sweep (1 sweep = $N$ flip attempts) and once the system reaches equilibrium we calculate the desired quantity. In order to determine whether the system has reached equilibrium or not we can look at the change in magnetization or internal energy with each sweep as both the quantities become stable on equilibrium.
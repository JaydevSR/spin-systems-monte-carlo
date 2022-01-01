#=
Exact solutions of 2D (square) Ising model by Onsager

( Taking k_B = 1 and J=1 )
=#

using SpecialFunctions

const critical_T = 2 / log1p(sqrt(2))

"""
    exact_internal_energy(T::Float64)

Calculate the exact value of internal energy per site at temperature T.
"""
function exact_internal_energy(T::Float64)
    β = 1/T
    k = 2tanh(2β) / cosh(2β)
    j = 2tanh(2β)^2 - 1
    K_k = SpecialFunctions.ellipk(k^2)
    return -coth(2β) * (1 + (2/π) * j * K_k)
end

"""
    exact_specific_heat(T::Float64)

Calculate the exact value of specific hear per site at temperature T.
"""
function exact_specific_heat(T::Float64)
    β = 1/T
    k = 2tanh(2β) / cosh(2β)
    j = 2tanh(2β)^2 - 1
    K_k = SpecialFunctions.ellipk(k^2)
    E_k = SpecialFunctions.ellipe(k^2)
    return β^2 * coth(2β)^2 * (2/π) * (((j - 1//2)^2 + 7//4) * K_k - 2E_k - (1 - j) * π / 2)
end
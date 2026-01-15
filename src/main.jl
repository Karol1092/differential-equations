include("integrate.jl")

f(x) = x^2
n = 3
x, w = gauss_legendre_nodes_weights(n)

I = sum(w .* f.(x))
println("Całka na [-1,1]: ", I)
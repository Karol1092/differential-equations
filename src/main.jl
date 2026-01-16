using Plots, LinearAlgebra
include("integrate.jl")

gr()

function k(x)
    return x <= 1.0 ? 0.5 : 1.0
end

function solve(n)
    L = 3.0
    h = L / n
    nodes = collect(0:h:3.0)

    K = zeros(n + 1, n + 1)
    F = zeros(n + 1)

    for e in 1:n
        i, j = e, e + 1
        xa, xb = nodes[i], nodes[j]

        e1 = -1/h
        e2 = 1/h

        K[i, i] += integrate(x -> k(x) * e1 * e1, xa, xb)
        K[i, j] += integrate(x -> k(x) * e2 * e1, xa, xb)
        K[j, i] += integrate(x -> k(x) * e1 * e2, xa, xb)
        K[j, j] += integrate(x -> k(x) * e2 * e2, xa, xb)

        F[i] += integrate(x -> -k(x) * 1.0 * e1, xa, xb)
        F[j] += integrate(x -> -k(x) * 1.0 * e2, xa, xb)
    end

    K[1, 1] += 0.5
    F[1] -= 0.5

    K[n + 1, :] .= 0.0
    K[n + 1, n + 1] = 1.0
    F[n + 1] = 0.0

    w = K \ F

    u = w .+ nodes
    
    # display(K)
    # display(F)
    # println(u)

    p = plot(nodes, u, label="u(x) - MES", marker=:circle, gui = true)
    display(p)
    println("Press enter to close the window...")
    readline()
end

solve(3)
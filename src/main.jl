using Plots, LinearAlgebra
include("integration.jl")

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

    de1 = -1/h
    de2 = 1/h

        for e in 1:n
            i, j = e, e + 1
            xa, xb = nodes[i], nodes[j]

            K[i, i] += integrate_split(x -> k(x) * de1 * de1, xa, xb)
            K[i, j] += integrate_split(x -> k(x) * de2 * de1, xa, xb)
            K[j, i] += integrate_split(x -> k(x) * de1 * de2, xa, xb)
            K[j, j] += integrate_split(x -> k(x) * de2 * de2, xa, xb)

            F[i] += integrate_split(x -> -k(x) * 1.0 * de1, xa, xb)
            F[j] += integrate_split(x -> -k(x) * 1.0 * de2, xa, xb)
        end

    K[1, 1] += 0.5
    F[1] -= 0.5

    K[n + 1, :] .= 0.0
    K[n + 1, n + 1] = 1.0
    F[n + 1] = 0.0

    w = K \ F

    u = w .+ nodes

    p = plot(nodes, u, label="u(x) - MES", marker=:circle, gui = true)
    display(p)
    println("Press enter to close the window...")
    readline()
end

print("n = ")
n = parse(Int, readline())
solve(n)
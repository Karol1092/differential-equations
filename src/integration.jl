function integrate(f, a, b)
    nodes = [-1 / sqrt(3), 1 / sqrt(3)]
    weights = [1.0, 1.0]
    
    # x_mapped = mx + k
    m = (b - a) / 2
    k = (a + b) / 2
    
    sum_val = 0.0
    for i in 1:2
        x_mapped = m * nodes[i] + k
        sum_val += weights[i] * f(x_mapped)
    end
    
    return sum_val * m
end

function integrate_split(f, ax, bx)
    if ax < 1 < bx
        return integrate(f, ax, 1) + integrate(f, 1, bx)
    else
        return integrate(f, ax, bx)
    end
end
function legendreP(n, x)
    if n == 0
        return 1.0
    elseif n == 1
        return x
    else
        P0, P1 = 1.0, x
        for k in 1:n-1
            Pk1 = ((2k+1)*x*P1 - k*P0)/(k+1)
            P0, P1 = P1, Pk1
        end
        return P1
    end
end

function legendrePprime(n, x)
    return n/(1 - x^2) * (legendreP(n-1, x) - x*legendreP(n, x))
end

function gauss_legendre_nodes_weights(n; tol=1e-14)
    x = zeros(n)
    w = zeros(n)
    
    for i in 1:n
        x0 = cos(pi*(i - 0.25)/(n + 0.5))
        xi = x0
        while true
            dx = legendreP(n, xi)/legendrePprime(n, xi)
            xi -= dx
            if abs(dx) < tol
                break
            end
        end
        x[i] = xi
        w[i] = 2 / ((1 - xi^2) * (legendrePprime(n, xi)^2))
    end
    
    return x, w
end


function integrate(f, a, b)
    n = 3
    x, w = gauss_legendre_nodes_weights(n)
    x_scaled = (b-a)/2 .* x .+ (a+b)/2
    w_scaled = (b-a)/2 .* w
    I2 = sum(w_scaled .* f.(x_scaled))
    return I2
end

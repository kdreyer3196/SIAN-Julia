include("IdentifiabilityODE.jl")

println("Setting up the problem")

R, (x1, x2, x3, x4, y1, y2, beta, c, d, k1, k2, mu1, mu2, q1, q2, s) = Nemo.PolynomialRing(Nemo.QQ, ["x1", "x2", "x3", "x4", "y1", "y2", "beta", "c", "d", "k1", "k2", "mu1", "mu2", "q1", "q2", "s"])

sigma_x = [
            [x1, -beta * x1 * x4 - d * x1 + s],
            [x2, beta * q1 * x1 * x4 - k1 * x2 - mu1 * x2],
            [x3, beta * q2 * x1 * x4 + k1 * x2 - mu2 * x3],
            [x4, -c * x4 + k2 * x3]
          ]

sigma_y = [
            [y1, x1],
            [y2, x4]
          ]

identifiability_ode(sigma_x, sigma_y,  [], get_parameters(sigma_x, sigma_y, []); p = 0.99, p_mod = 2^29 - 3, nthrds = 64)


include("../src/symfunc.jl")

# Parameters; see Osipenko (2006)
r = 2.0
a = -0.9
b = 0.9
C1 = 0.4
C3 = 6.0

# Symbolic computations
x = Sym("x")
y = Sym("y")
tau = C1 - C3 / (1.0 + x^2 + y^2)
fx = r + a * (x * cos(tau) - y * sin(tau))
fy = b * (x * sin(tau) + y * cos(tau))
symfunc_2_juliafunc((fx, fy), (x, y), base_name = "Ikeda")

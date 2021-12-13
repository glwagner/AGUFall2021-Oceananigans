using Oceananigans, Oceananigans.Units

# Grid
Nx = Ny = 60
Lx = Ly = 1200kilometers

grid = RectilinearGrid(size = (Nx, Ny),
                       x = (-Lx/2, Ly/2),
                       y = (-Ly/2, Ly/2),
                       topology = (Bounded, Bounded, Flat))

# Coriolis and viscosity
coriolis = BetaPlane(f₀=1e-4, β=1e-11)
closure = IsotropicDiffusivity(ν=4e2)

# Top boundary condition
wind_stress(x, y, t) = -1e-4 * cos(π * x / 600kilometers)
u_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(wind_stress))

# Model
model = NonhydrostaticModel(; grid, coriolis, closure,
                            boundary_conditions = (; u=u_bcs))
                           
simulation = Simulation(model, Δt = 3600, stop_time = 1years)

run!(simulation)


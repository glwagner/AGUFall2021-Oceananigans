using Oceananigans, Oceananigans.Units, GLMakie

# Specify domain and mesh / grid
grid = RectilinearGrid(size = (60, 60, 1),
                       x = (0, 1200kilometers),
                       y = (0, 1200kilometers),
                       z = (-5000, 0),
                       halo = (3, 3, 3),
                       topology = (Bounded, Bounded, Bounded))

# Specify coriolis, viscosity / turbulence closure, and boundary condition
coriolis = BetaPlane(f₀=1e-4, β=1e-11)
closure = IsotropicDiffusivity(ν=400)

wind_stress(x, y, t) = - 1e-4 * cos(π * x / 600kilometers)
u_bcs = FieldBoundaryConditions(; top=FluxBoundaryCondition(wind_stress))

model = HydrostaticFreeSurfaceModel(; grid, coriolis, closure,
                                      momentum_advection = UpwindBiasedFifthOrder(),
                                      boundary_conditions = (; u=u_bcs))
                          
# Time-stepping
simulation = Simulation(model, Δt=20minutes, stop_time=3years)

progress(sim) = @info "Iter: $(iteration(sim)), time: $(prettytime(sim))"
simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

run!(simulation)

fig, ax, pl = heatmap(interior(model.velocities.u)[:, :, 1])
display(fig)

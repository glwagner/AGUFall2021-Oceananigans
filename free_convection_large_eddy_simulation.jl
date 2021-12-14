using Oceananigans, Oceananigans.Units, GLMakie

grid = RectilinearGrid(GPU(), size = (256, 256, 256), halo = (3, 3, 3),
                       x = (0, 256), y = (0, 256), z = (-128, 0),
                       topology = (Periodic, Periodic, Bounded))

# Specify coriolis, viscosity / turbulence closure, and boundary condition
coriolis = FPlane(f=1e-4)
closure = AnisotropicMinimumDissipation()
boundary_conditions = (; b = FieldBoundaryConditions(top=FluxBoundaryCondition(1e-8)))

model = NonhydrostaticModel(; grid, coriolis, closure, boundary_conditions,
                              advection = WENO5(), tracers = :b,
                              buoyancy = BuoyancyTracer())

set!(model, b = (x, y, z) -> 1e-6 * z, w = (x, y, z) -> 1e-4 * randn())
                          
# Time-stepping
simulation = Simulation(model, Δt=10, stop_time=12hours)

progress(sim) = @info "Iter: $(iteration(sim)), time: $(prettytime(sim)), " *
                      "Δt: $(prettytime(sim.Δt)), wall time: $(prettytime(sim.run_wall_time))"
simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

wizard = TimeStepWizard(cfl=0.2, max_change=1.1)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(100))

times, slices = [], []
function slice_w(sim)
    push!(times, sim.model.clock.time)
    push!(slices, Array(interior(sim.model.velocities.w))[:, 1, :])
    return nothing
end
simulation.callbacks[:slice_catcher] = Callback(slice_w, TimeInterval(5minute))

run!(simulation)

# Visualize vertical velocity
wmax = maximum(abs, model.velocities.w)
x, y, z = nodes(model.velocities.w)
n = Node(1)
fig = Figure(resolution=(1200, 900))
title = @lift "Vertical velocity at t = " * prettytime(times[$n])
ax = Axis(fig[1, 1], title=title, xlabel="x (m)", ylabel="z (m)")
hm = heatmap!(ax, x, z, @lift(slices[$n]), colormap=:balance, limits=(-wmax, wmax))
Colorbar(fig[1, 2], hm)

record(fig, "free_convection.mp4", 1:length(times), framerate=12) do nn
    n[] = nn
end

export SimulationInstance

struct SimulationInstance
    setup_config::Function
    setup_benchmark::Function
    create_mesh::Function
    create_simulate::Function
    run_simulation::Function
end

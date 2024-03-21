using Random, BrkgaMpIpr

include("osvrp_instance.jl")
include("osvrp_decoder.jl")
include("osvrp_greedy.jl")

########################################
# Read arguments and Instance
########################################

# p_file = "Instâncias\\GP03-01.txt"
# instance = OSVRPInstance(p_file)
# cost = BRKGA(p_file, 2, 2, 0.30, 0.15, 0.55, 10, 1000, 10)
#α=10;β=10;pe=0.30;pm=0.22;ρ=0.55;L=10;R=1000;λ=10

function BRKGA(p_file::String, α::Int64, β::Int64, pe::Float64, pm::Float64, ρ::Float64, L::Int64, R::Int64, λ::Int64)

    instance = OSVRPInstance(p_file)

    Size = length(instance.ρ) + instance.n + instance.n

    seed = rand(1:1000)
    configuration_file = "config.conf"
    num_generations = β*Size

    brkga_params, control_params = load_configuration(configuration_file)

    brkga_params.population_size = α*Size
    brkga_params.elite_percentage = pe
    brkga_params.mutants_percentage = pm

    ########################################
    # Build the BRKGA data structures and initialize
    ########################################

    brkga_data = build_brkga(instance, osvrp_decode_ND!, MINIMIZE, seed, Size,
        brkga_params, true
    )

    # Set Bias probability
    set_bias_custom_function!(brkga_data, x -> x == 1 ? ρ : 1.0 - ρ)

    # Warm-start
    seq_spt, route_spt = NEH_VRP(instance, false)
    key_spt = osvrp_encode_ND!(seq_spt, route_spt, instance)

    seq_lpt, route_lpt = NEH_VRP(instance, true)
    key_lpt = osvrp_encode_ND!(seq_lpt, route_lpt, instance)

    set_initial_population!(brkga_data, [key_spt, key_lpt])

    # NOTE: don't forget to initialize the algorithm.
    initialize!(brkga_data)

    # Warmup - Algorithm BRKGA
    bogus_data = deepcopy(brkga_data)
    evolve!(bogus_data, 2)

    get_best_fitness(brkga_data)
    get_best_chromosome(brkga_data)
    bogus_data = nothing

    bestsol = minimum([osvrp_decode_ND!(key_spt, instance, true),
                    osvrp_decode_ND!(key_spt, instance, true)])

    stableit = 0

    ########################################
    # Find good solutions / evolve
    ########################################

    t0 = time_ns()
    for gen in 1:num_generations
        #global bestsol
        evolve!(brkga_data) #Evolve P one generation!

        actualsol = get_best_chromosome(brkga_data)

        if mod(gen, L) == 0 # Realizar a busca local
            seq, route = osvrp_decodeND(actualsol, instance)
            seq, route = Insertion_Local_Search(seq, route, instance)

            key_ls = osvrp_encode_ND!(seq, route, instance)

            inject_chromosome!(brkga_data, key_ls, 1, 1)
        end

        cursol = get_best_fitness(brkga_data)

        if cursol == bestsol
            stableit = stableit + 1
            if stableit ≥ R
                actualsol = get_best_chromosome(brkga_data)
                reset!(brkga_data)
                inject_chromosome!(brkga_data, actualsol, 1, 1)
                inject_chromosome!(brkga_data, key_spt, 1, 2)
                inject_chromosome!(brkga_data, key_lpt, 1, 3)
                stableit = 0
            end

        else
            stableit = 0
        end

        if cursol < bestsol
            bestsol = cursol
        end

        if stableit == 10000
            break
        end

        if (time_ns() - t0)*1e-9 >= 1800
            break
        end
    end

    cost = bestsol
    println("Instância $(rsplit(p_file, "\\")[2]) -- Custo: $cost")
    return bestsol
end

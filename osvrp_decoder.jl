function completion_time(oper_sequence::Array{Int64}, instance::OSVRPInstance)
    M, J = zeros(Int, instance.m), zeros(Int, instance.n)

    for idx in eachindex(oper_sequence)
        operation = oper_sequence[idx]

        machine, job = instance.O[operation]
        process_time = instance.ρ[machine, job]

        M[machine] = max(M[machine], J[job]) + process_time
        J[job] = M[machine]
    end

    return J
end

function completion_time_m_AS(oper_sequence::Array{Int64},
                                instance::OSVRPInstance)

    M, J = zeros(Int, instance.m), zeros(Int, instance.n)
    C = copy(oper_sequence)
    O = instance.O

    while length(C) > 0

        st = [max(M[O[operation][1]], J[O[operation][2]]) for operation in C]
        yt = [max(M[O[operation][1]], J[O[operation][2]]) + instance.ρ[O[operation][1], O[operation][2]] for operation in C]

        y = minimum(st)
        z = minimum(yt)

        maxStartTime = floor(z - 0.5*(z-y))
        Gt = [operation for (index, operation) in enumerate(C) if st[index]<=maxStartTime]
        early_operation = Gt[1]

        filter!(x->x≠early_operation, C)

        machine, job = O[early_operation]
        process_time = instance.ρ[machine, job]

        M[machine] = max(M[machine], J[job]) + process_time
        J[job] = M[machine]
    end

    return J
end

function  time_delivery(route::Array{Int64}, J::Array{Int64}, instance::OSVRPInstance)
    n = instance.n - 1

    R = J[route[1]] + instance.δ[1, route[1]+1]
    A = instance.σ[route[1]]
    t_exit = J[route[1]]

    for i in 1:n
        idy = route[i]
        global idx = route[i+1]
        if (A + instance.σ[idx]) <= instance.ϕ && R >= J[idx] && t_exit >= J[idx]
            R = R + instance.δ[idy+1, idx+1]
            A += instance.σ[idx]
        else # Voltando para o Depot
            R = R + instance.δ[idy+1, 1]
            R = max(R, J[idx])
            t_exit = R
            A = 0
            R = R + instance.δ[1, idx+1]
            A += instance.σ[idx]
        end
    end

    R += instance.δ[idx+1, 1]

    return R
end

function  routing(job_seq::Any, route_idx::Any, J::Array{Int64},
                                    instance::OSVRPInstance)
    m = maximum(route_idx)#Number of routes
    route = [[job for job in job_seq if route_idx[job]==route] for route in 1:m]
    route = [sub_rot for sub_rot in route if length(sub_rot) > 0]

    C_max = (sub_route) -> maximum([J[i] for i in sub_route])
    sort!(route, by=C_max)

    makespan = 0
    for sub_route in route
        for (idx, job) in enumerate(sub_route)

            if idx==1
                makespan = maximum((makespan, C_max(sub_route)))
                makespan += instance.δ[1, job+1]
            end
            if idx +1 <= length(sub_route)
                makespan += instance.δ[job+1, sub_route[idx+1]+1]
            end
            if idx == length(sub_route)
                makespan += instance.δ[job+1, 1]
            end
        end
    end

    # Penalizando a solução se existir rotas inviáveis
    for sub_route in route
        sum_sub = sum(Float64[instance.σ[j] for j in sub_route])
        if sum_sub > instance.ϕ
            makespan += (sum_sub - instance.ϕ)*instance.ϕ
        end
    end

    return makespan
end

function osvrp_decode!(chromosome::Array{Float64}, instance::AbstractInstance,
                     rewrite::Bool)::Int64

    quant_operations = length(instance.ρ)

    oper_sequence = chromosome[1:quant_operations]
    oper_sequence = sortperm(oper_sequence)# Sequência de processamento das operações

    route = chromosome[quant_operations+1:end]
    route = sortperm(route)# Sequência de entrega das ordens para os clientes

    J = completion_time(oper_sequence, instance)# Pegando os tempos de términos das tarefas
    makespan = time_delivery(route, J, instance)# Roteirizando as tarefas

    return makespan
end

function osvrp_decode_ND!(chromosome::Array{Float64}, instance::AbstractInstance,
                     rewrite::Bool)::Int64

    m = instance.n#ceil(Int, sum(instance.σ)/instance.ϕ)+2

    quant_operations = length(instance.ρ)

    oper_sequence = chromosome[1:quant_operations]
    oper_sequence = sortperm(oper_sequence)# Sequência de processamento das operações

    job_seq = chromosome[quant_operations+1:end-instance.n]
    job_seq = sortperm(job_seq)# Sequência de entrega das ordens para os clientes

    route_idx = chromosome[end-instance.n+1:end]
    route_idx = [maximum([1.0, ceil(Int, m*idx)]) for idx in route_idx]

    J = completion_time(oper_sequence, instance)# Pegando os tempos de términos das tarefas
    makespan = routing(job_seq, route_idx, J, instance)# Roteirizando as tarefas

    return makespan
end

function osvrp_encode_ND!(seq::Array{Int64}, route::Any, instance::AbstractInstance)
    route = [rot for rot in route if length(rot) > 0]
    m = length(route)

    route_seq = [job for sub_route in route for job in sub_route if length(sub_route) > 0]
    route_idx = zeros(instance.n)

    for job in 1:instance.n
        for (idx, sub_route) in enumerate(route)
            if job in sub_route
                route_idx[job] = idx
            end
        end
    end

    Size = length(instance.ρ) + instance.n + instance.n
    keys = sort(rand(Size))
    initial_chromosome = zeros(Size)

    for i in 1:length(instance.ρ)
        initial_chromosome[seq[i]] = keys[i]
    end

    for i in 1:instance.n
        initial_chromosome[route_seq[i]+length(instance.ρ)] = keys[i+length(instance.ρ)]
    end

    for job in 1:instance.n
        initial_chromosome[job+length(instance.ρ)+instance.n] = route_idx[job]/m
    end
    return initial_chromosome
end

function osvrp_decodeND(chromosome::Array{Float64}, instance::AbstractInstance)

    m = instance.n#ceil(Int, sum(instance.σ)/instance.ϕ)+2

    quant_operations = length(instance.ρ)

    oper_sequence = chromosome[1:quant_operations]
    oper_sequence = sortperm(oper_sequence)# Sequência de processamento das operações

    job_seq = chromosome[quant_operations+1:end-instance.n]
    job_seq = sortperm(job_seq)# Sequência de entrega das ordens para os clientes

    route_idx = chromosome[end-instance.n+1:end]
    route_idx = [maximum([1.0, ceil(Int, m*idx)]) for idx in route_idx]

    route = [[] for rot in 1:m]
    for r_idx in 1:m
        for job in job_seq
            if route_idx[job] == r_idx
                push!(route[r_idx], job)
            end
        end
    end
    return oper_sequence, route
end
#include("osvrp_instance.jl")
#p_file = "Instâncias\\GP03-01.txt"
#instance = OSVRPInstance(p_file)
# seq = [27, 29, 14, 24, 7, 11, 13, 33, 1, 16, 20, 21, 19, 30, 32, 34, 23, 2, 9, 10, 31, 36, 5, 6, 35, 25, 15, 3, 18, 17, 8, 26, 12, 22, 4, 28]
# route =  [[], [3, 4, 2, 5, 1], [6]]
#
# key = osvrp_encode_ND!(seq, route, instance)
# key_spt = [0.7110472378313153, 0.3544932419079556, 0.022428049107474424, 0.28377881692070606, 0.06699028837167353, 0.09601634406393389, 0.6899231779594881, 0.24534510241402585, 0.6914031822771953, 0.20647263466623222, 0.4893088829873495,
# 0.5804155276082061, 0.4692464464553492, 0.6252803944746839, 0.30271089250984407, 0.48242861926877234, 0.014996949947458038, 0.045171873344632685, 0.1504877315169033, 0.3519083484624521, 0.6312873472866574, 0.21837945146943527,
#0.6731569333686342, 0.49387041475001836, 0.339920480708062, 0.4478692191089282, 0.3008831265491996, 0.35724858326009845, 0.7272439135608377, 0.012034292356487564, 0.16452637723622843, 0.6139219799701565, 0.6004751387547262, 0.6137787310424501, 0.1964133064841742, 0.7282693389109918, 0.8480303282671395, 0.7652376198878714, 0.850232148987919, 0.8262866640454638, 0.7660585561878679, 0.7916682934114725, 1.0, 0.5, 1.0, 1.0, 0.5, 0.5]

#key_paper = [0.01, 0.46, 0.24, 0.65, 0.37, 0.22, 0.43, 0.15, 0.68, 0.79, 0.77, 0.82, 0.65, 0.33, 0.56]
#s, r = osvrp_decodeND(key_paper, instance)

# seq, route = NEH_VRP(instance, false)
#
# key = osvrp_encode_ND!(seq, route, instance)
#
# make = osvrp_decode_ND!(key, instance, true)

# p_file = "Instâncias\\GP03-01.txt"
# instance = OSVRPInstance(p_file)
#
#key = [0.4, 0.2, 0.6, 0.1, 0.8, 0.3, 0.9, 0.5, 0.7, 0.1, 0.3, 0.2, 0.40, 0.60, 0.4, 0.4, 0.2, 0.4]
#objective = osvrp_decode_ND!(key, instance, true)
#s, r = osvrp_decodeND(key, instance)

# p_file = "Instâncias\\GP05-03.txt"
# instance = OSVRPInstance(p_file)
#
# J = [1023, 1257, 1265, 1024, 1023]
# route = [1, 4, 5, 2, 3]
#
# n = instance.n - 1
#
# R = J[route[1]] + instance.δ[1, route[1]+1]
# A = instance.σ[route[1]]
# t_exit = J[route[1]]
#
# for i in 1:n
#     global A, R, t_exit
#     idy = route[i]
#     global idx = route[i+1]
#     if (A + instance.σ[idx]) <= instance.ϕ && R >= J[idx] && t_exit >= J[idx]
#         R = R + instance.δ[idy+1, idx+1]
#         A += instance.σ[idx]
#     else # Voltando para o Depot
#         R = R + instance.δ[idy+1, 1]
#         A = 0
#         R = max(R, J[idx]) + instance.δ[1, idx+1]
#         A += instance.σ[idx]
#         t_exit = J[idx]
#     end
# end
#
# R += instance.δ[idx+1, 1]
#
# println(R)

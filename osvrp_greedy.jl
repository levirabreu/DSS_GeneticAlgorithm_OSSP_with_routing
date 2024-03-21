function  routing_NEH(route::Any, J::Array{Int64},
                                    instance::OSVRPInstance)
    route = [sub_rot for sub_rot in route if length(sub_rot) > 0]

    if length(route) == 0 # Não possui rota, a solução é o makespan atual
        return maximum(J)
    end

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

    return makespan
end

function completion_time_seq(oper_sequence::Array{Int64}, instance::OSVRPInstance)
    M, J, JM = zeros(Int, instance.m), zeros(Int, instance.n), zeros(Int, instance.n)

    for idx in eachindex(oper_sequence)
        operation = oper_sequence[idx]

        machine, job = instance.O[operation]
        process_time = instance.ρ[machine, job]

        M[machine] = max(M[machine], J[job]) + process_time
        J[job] = M[machine]
        JM[job] += 1
    end

    return J, JM
end

function completion_time_m_AS_seq(oper_sequence::Array{Int64},
                                instance::OSVRPInstance)

    M, J, JM = zeros(Int, instance.m), zeros(Int, instance.n), zeros(Int, instance.n)
    C = copy(oper_sequence)
    O = instance.O
    S = []

    while length(C) > 0

        st = [max(M[O[operation][1]], J[O[operation][2]]) for operation in C]

        yt = [max(M[O[operation][1]], J[O[operation][2]]) + instance.ρ[O[operation][1], O[operation][2]] for operation in C]

        y = minimum(st)
        z = minimum(yt)

        maxStartTime = floor(z - 0.5*(z-y))
        Gt = [operation for (index, operation) in enumerate(C) if st[index]<=maxStartTime]
        early_operation = Gt[1]

        filter!(x->x≠early_operation, C)
        push!(S, early_operation)


        machine, job = O[early_operation]
        process_time = instance.ρ[machine, job]

        M[machine] = max(M[machine], J[job]) + process_time
        J[job] = M[machine]
        JM[job] += 1
    end

    return J, JM, S
end

function saving(k::Int64, i::Int64, j::Int64, instance::OSVRPInstance)::Int64
    Δ_length = instance.δ[i+1, k+1] + instance.δ[k+1, j+1] - instance.δ[i+1, j+1]
    return Δ_length
end

function NEH_VRP(instance::OSVRPInstance, lpt::Bool=false)

    Size = instance.m*instance.n
    p_idx = x -> instance.ρ[instance.O[x][1], instance.O[x][2]]
    seq_initial = sort(collect(1:Size), by=p_idx, rev=lpt)
    seq_initial = completion_time_m_AS_seq(seq_initial, instance)[3]

    s = Int64[]

    Δ_length = 0
    best_pos_perm = 0
    best_sub_route = 0
    best_pos_route = 0
    best_pos_r = 0
    JM = []
    route = [[] for n in 1:instance.n]

    while length(s) < Size
        operation = seq_initial[1]
        job = instance.O[operation][2]
        best_make = Inf
        for idx in 1:length(s)+1
            ss = insert!(copy(s), idx, operation)
            J, JM = completion_time_seq(ss, instance)

            if JM[job] == instance.m # Se o job já tiver passado em todas as maq.
                for (id_r, sub_route) in enumerate(route)
                    if instance.σ[job] + sum(Float64[instance.σ[j] for j in sub_route]) > instance.ϕ
                        continue # A rota não possui capacidade para o job
                    end

                    # Encontrando a melhor posição para colocar o job na subrota
                    best_length = Inf
                    best_pos_r = 0
                    for pos_r in 1:length(sub_route)+1
                        if pos_r == 1
                            Δ_length =  saving(job, 0, 0, instance)
                        elseif pos_r == length(sub_route)+1
                            Δ_length = saving(job, sub_route[pos_r-1], 0, instance)
                        else
                            Δ_length = saving(job, sub_route[pos_r-1], sub_route[pos_r], instance)
                        end

                        if Δ_length < best_length
                            best_length = Δ_length
                            best_pos_r = pos_r
                        end
                    end

                    rr = copy(route)
                    rr[id_r] = insert!(copy(sub_route), best_pos_r, job)
                    make_ss = routing_NEH(rr, J, instance)

                    if make_ss < best_make
                        best_make = make_ss
                        best_pos_perm = idx
                        best_sub_route = id_r
                        best_pos_route = best_pos_r
                    end
                end
            else
                make_ss = routing_NEH(route, J, instance)

                if make_ss < best_make
                    best_make = make_ss
                    best_pos_perm = idx
                end
            end
        end

        if JM[job] == instance.m # Job já passou em todas as máquinas
            insert!(s, best_pos_perm, operation)
            insert!(route[best_sub_route], best_pos_route, job)
            popfirst!(seq_initial)
        else
            insert!(s, best_pos_perm, operation)
            popfirst!(seq_initial)
        end
    end
    return s, route
end

function Insertion_Local_Search(seq::Array{Int64}, route::Any, instance::OSVRPInstance)
    n_operations = instance.n*instance.m
    improve = true

    seq_best = copy(seq)
    route_best = copy(route)
    best_make = routing_NEH(route_best, completion_time_seq(seq_best, instance)[1], instance)

    best_ss = 0
    best_rr = 0

    while improve == true
        improve = false
        oper_random = randperm(n_operations)
        for oper in oper_random
            seq_aux = filter!(x->x≠oper, copy(seq_best))
            job = instance.O[oper][2]
            route_aux = [[k for k in sub_route if k!= job] for sub_route in route_best]
            improve_aux = false
            for idx in 1:length(seq_aux)+1
                ss = insert!(copy(seq_aux), idx, oper)
                J, JM = completion_time_seq(ss, instance)

                # Testar o job em todas as posições pegar a melhor!

                for (r_idx, sub_route) in enumerate(route_aux)
                    for pos_r in 1:length(sub_route)+1
                        rr = copy(route_aux)
                        rr[r_idx] = insert!(copy(sub_route), pos_r, job)

                        make = routing_NEH(rr, J, instance)

                        if make < best_make
                            best_make = make
                            improve = true
                            improve_aux = true
                            best_rr = copy(rr)
                            best_ss = copy(ss)
                        end
                    end
                end
            end

            if improve == true && improve_aux == true
                seq_best = copy(best_ss)
                route_best = copy(best_rr)
                best_make = routing_NEH(route_best, completion_time_seq(seq_best, instance)[1], instance)
            end
        end
    end
    return seq_best, route_best
end

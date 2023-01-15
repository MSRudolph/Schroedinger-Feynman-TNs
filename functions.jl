using ITensors
include("datatypes.jl")
include("gates.jl")
include("circuits.jl")

function get_indices(n_qubits::Int)
    # return [Index(2, "qubit"*string(ii)) for ii in 1:n_qubits]
    return siteinds("S=1/2", n_qubits)
end


function get_initial_wavefunction(indices)
    n_qubits = length(indices)
    
    wf = zeros(2^n_qubits)
    wf[1] = 1
    wf = reshape(wf, [2 for _ in 1:n_qubits]...)
    return Wavefunction(itensor(wf, indices), indices)

end

function get_initial_wavefunction(indices::Index)
    return get_initial_wavefunction([indices])
end

function get_initial_clusterwavefunction(indices, splits)
    return ClusterWavefunction(
        [get_initial_wavefunction(indices[qinds]) for qinds in splits], 
        indices
        )
end

function get_initial_mps(indices)
    return productMPS(indices, "Up")
end


################    apply    #################################

function apply(state, gates::Vector{Gate}; kwargs...)
    for gate in gates
        state = apply(state, gate; kwargs...)
    end
    return state
end

function apply(states::Vector, gate::Gate; kwargs...)
    all_states::Vector{eltype(states)} = []
    for state in states
        state = apply(state, gate; kwargs...)
        all_states = add_to_list!(all_states, state)
    end
    return all_states
end

function apply(states, circuits::Vector{Vector{Gate}}; kwargs...)
    all_states::Vector{typeof(states)} = []
    for circuit in circuits
        state = apply(states, circuit; kwargs...)
        add_to_list!(all_states, state)
    end
    return all_states
end


function apply(states::Vector, gate::Gate, maxreplica; kwargs...)
    start_projecting::Bool = length(states) >= maxreplica
    states = apply(states, gate; start_projecting=start_projecting, kwargs...)
    if start_projecting
        states = states[1:floor(Int, maxreplica)]
    end
    return states
end


### Exact ###

function apply(wavefunction::Wavefunction, gate::Gate; kwargs...)
    return Wavefunction(noprime!(wavefunction.itensor * gate.itensor), wavefunction.indices)
end


### Cluster ###

function apply(state::ClusterWavefunction, gates::Vector{Gate}, maxreplica; kwargs...)
    for gate in gates
        state = apply(state, gate, maxreplica; kwargs...)
    end
    return state
end

function apply(state::ClusterWavefunction, gate::Gate; start_projecting=false, kwargs...)
    active_wf_inds = find_active_indices(state, gate)
    ##
    if allequal(active_wf_inds)
        # contract the gate with one wavefunction
        state = copy(state)
        state.wavefunctions[active_wf_inds[1]] = apply(state.wavefunctions[active_wf_inds[1]], gate; kwargs...)
        return state
    else
        # this must be a two-qubit gate on two different wavefunctions
        # decompose the gate, replicate the ClusterWavefunction, and apply the sets of gates to the corresponding wavefunctions
        gates = decompose(gate)
        if start_projecting
            gates = gates[rand(1:2)]
        end
        states = apply(state, gates)
        
        return states 
    end
end

function apply(state::ClusterWavefunction, gate::Gate, maxreplica; kwargs...)
    return apply(state, gate; kwargs...)
end


### MPS ###

function apply(state::MPS, gate::Gate; cutoff=1e-10, maxdim=10_000, maxdist=3)

    if (length(gate.itensor.tensor.inds) >= 4) && (maxdist<length(state)) # if multi-qubit gates
        
        relevant_ids = gate.itensor.tensor.inds[1:2]

        sites_found::Vector{Int} = [findall(x->x==ind, siteinds(state))[1] for ind in relevant_ids]

        if diff(sites_found)[1] > maxdist
            gates = decompose(gate)
            return apply(state, gates; cutoff=cutoff, maxdim=10_000, maxdist=maxdist)
        end
    end
    return ITensors.apply(gate.itensor, state; cutoff=cutoff, maxdim=maxdim)
end



###############    Utilities    #############################

function as_statevector(states)
    wf = as_wavefunction(states)
    return as_statevector(wf)
end

function as_statevector(wf::Wavefunction)
    return reshape(array(wf.itensor, wf.indices), 2^length(wf.indices))
end

function as_statevector(wfs::Vector{Wavefunction})
    return sum(as_statevector.(wfs))
end

function as_wavefunction(states)
    return as_wavefunction.(states)
end

function as_wavefunction(cluster_wf::ClusterWavefunction)
    return Wavefunction(itensor(array(contract([wf.itensor for wf in cluster_wf.wavefunctions]), cluster_wf.indices), cluster_wf.indices), cluster_wf.indices)
end

function as_wavefunction(wf::Wavefunction)
    return finalize(wf)
end

function as_wavefunction(mps::MPS)
    return Wavefunction(contract(mps), siteinds(mps))
end


function fidelity(state1, state2)
    return abs(as_statevector(state1)' * as_statevector(state2))
end




###############    Helpers    ###########################

using StatsBase
function remove!(states, fidelity)
    if (length(states) > 1)
        rm_inds = findall(x->x==1, StatsBase.wsample([0, 1], [fidelity, 1-fidelity], length(states); replace=true))
        # @show rm_inds
        if 0 < length(rm_inds) < length(states)
            # println("Removing ", length(rm_inds))
            deleteat!(states, rm_inds)
        elseif length(rm_inds) == length(states)
            # println("Removing ", length(rm_inds)-1)
            deleteat!(states, rm_inds[1:end-1])
        end
    end
end


function add_to_list!(list, obj_to_add::Vector)
    return append!(list, obj_to_add)
end

function add_to_list!(list, obj_to_add)
    return push!(list, obj_to_add)
end

function finalize(wf::Wavefunction)
    return Wavefunction(itensor(array(wf.itensor, wf.indices), wf.indices), wf.indices)
end


function find_active_indices(state::ClusterWavefunction, gate)
    all_gate_inds = gate.itensor.tensor.inds
    gate_inds = all_gate_inds[1:length(all_gate_inds)/2]
    ##
    active_wf_inds = zeros(Int, length(gate_inds))
    for (ii, ind) in enumerate(gate_inds)
        for (kk, wf) in enumerate(state.wavefunctions)
            if ind in wf.indices
                active_wf_inds[ii] = kk
                continue
            end
        end
    end
    return active_wf_inds
end


function find_active_indices(splits, gate::Gate, indices)
    all_gate_inds = gate.itensor.tensor.inds
    gate_inds = all_gate_inds[1:length(all_gate_inds)/2]
    ##
    active_wf_inds = zeros(Int, length(gate_inds))
    for (ii, ind) in enumerate(gate_inds)
        ind = findall(x -> x==ind, indices)[1]
        for (kk, sp) in enumerate(splits)
            if ind in sp 
                active_wf_inds[ii] = kk
                continue
            end
        end
    end
    return active_wf_inds
end

function find_active_indices(splits, pair)
    ##
    active_wf_inds = zeros(Int, length(pair))
    for (ii, ind) in enumerate(pair)
        for (kk, sp) in enumerate(splits)
            if ind in sp 
                active_wf_inds[ii] = kk
                continue
            end
        end
    end
    return active_wf_inds
end
################ estimate memory ##########################

function memory_req(n_qubits, max_qubits_per_wf, depth; fidelity=1.0)
    # assumes linear topolgy
    n_wfs = ceil(n_qubits/max_qubits_per_wf)
    n_qubits_per_wf = partition_n(n_qubits, max_qubits_per_wf)
    n_cross_gates = n_wfs-1
    return memory_formula(n_cross_gates, n_qubits_per_wf, depth; fidelity=fidelity)
end

function memory_req(n_qubits, max_qubits_per_wf, depth, topology; fidelity=1.0)
    # takes any topolgy
    n_qubits_per_wf = partition_n(n_qubits, max_qubits_per_wf)
    # @show n_qubits_per_wf
    splits = get_splits_from_partition(n_qubits_per_wf)
    # @show topology, splits
    n_cross_gates = calculate_n_of_cross_gates(topology, splits)
    # @show n_cross_gates
    return memory_formula(n_cross_gates, n_qubits_per_wf, depth; fidelity=fidelity)
end

MEM_LIMIT::Float64 = 1e100

function memory_formula(n_cross_gates, n_qubits_per_wf, depth; fidelity=1.0)
    if maximum(n_qubits_per_wf) > 22 || (depth * n_cross_gates) > 62
        # @show n_qubits_per_wf, depth, n_cross_gates
        return MEM_LIMIT
    end
    # @show depth, n_cross_gates, 2^(depth * n_cross_gates)
    return ceil(fidelity * 2^(depth * n_cross_gates) * sum([2^nq for nq in n_qubits_per_wf])) 
end

function find_optimal_split(n_qubits, depth; fidelity=1.0)
    # assumes linear topolgy
    return _find_optimal_split(n_qubits, depth, get_line_topology(n_qubits); fidelity=fidelity)
end

function find_optimal_split(n_qubits, depth, topology; fidelity=1.0)
    # assumes linear topolgy
    return _find_optimal_split(n_qubits, depth, topology; fidelity=fidelity)
end

function _find_optimal_split(n_qubits, depth, topology; fidelity=1.0)
    # assumes linear topolgy
    steps = push!(collect(1:Int(ceil(n_qubits/2))), n_qubits)
    memory_estimates = [memory_req(n_qubits, m, depth, topology; fidelity=fidelity) for m in steps]
    if allequal([MEM_LIMIT, memory_estimates...])
        return 1:n_qubits
    end
    opt_n_per_wf = steps[argmin(memory_estimates)]
    splits = get_splits_from_partition(partition_n(n_qubits, opt_n_per_wf))
    return splits
end

function get_splits_from_partition(n_qubits_per_wf)
    n_qubits_per_wf = append!([0], n_qubits_per_wf)
    splits = [(sum(n_qubits_per_wf[1:ii])+1):sum(n_qubits_per_wf[1:ii+1]) for ii in 1:(length(n_qubits_per_wf)-1)]
    return splits
end

function partition_n(n_qubits, n_per_wf)
    rounded_n_per_wf = Int(ceil(n_per_wf))
    # @show n_qubits, n_per_wf, rounded_n_per_wf, Int(ceil(n_qubits/n_per_wf))
    nqubits_per_wf = [(rounded_n_per_wf*ii <= n_qubits) ? rounded_n_per_wf : (n_qubits % rounded_n_per_wf) for ii in 1:Int(ceil(n_qubits/rounded_n_per_wf))]
    if nqubits_per_wf[end] == 1
        pop!(nqubits_per_wf)
        nqubits_per_wf[end] += 1
    end
    return nqubits_per_wf
end


########### Circuits #######################

function expected_num_paths(circuit, splits, indices)
    n_xgates = calculate_n_of_cross_gates(circuit, splits, indices)
    return 2^n_xgates
end

function calculate_n_of_cross_gates(circuit, splits, indices)
    n_xgates = 0
    for gate in circuit
        inds = find_active_indices(splits, gate, indices)
        if !allequal(inds)
            n_xgates += 1
        end
    end
    return n_xgates
end

function calculate_n_of_cross_gates(topology, splits)
    n_xgates = 0
    for pair in topology
        inds = find_active_indices(splits, pair)
        if !allequal(inds)
            n_xgates += 1
        end
    end
    return n_xgates
end

function decompose(gates::Vector{Gate}, cluster_state)
    # create a list of circuit
    # loop over gates in circuit
        # add inside-gates to all list of circuits
        # decompose cross-gates
        # double list of circuits, add the gate halves to one of two copies
    
    all_circuits::Vector{Vector{Gate}} = [[]]
    for (ii, gate) in enumerate(gates)
        active_inds = find_active_indices(cluster_state, gate)
        if allequal(active_inds)
            add_to_circuit!(all_circuits, gate)
        else
            gates = decompose(gate)
            c1 = deepcopy(all_circuits)
            add_to_circuit!(c1, gates[1])
            
            c2 = all_circuits
            add_to_circuit!(c2, gates[2])
            all_circuits = [c1; c2]
        end
    end
    return all_circuits
end

function decompose(gates::Vector{Gate}, cluster_state, fidelity)
    all_circuits = decompose(gates, cluster_state)
    return StatsBase.sample(all_circuits, round(Int, fidelity*length(all_circuits)); replace=false)
end


function add_to_circuit!(circuit::Vector{Gate}, gate)
    return add_to_list!(circuit, gate)
end

function add_to_circuit!(circuits, gate)
    for (ii, circuit) in enumerate(circuits)
        circuits[ii] = add_to_list!(circuit, gate)
    end
    return circuits
end
 

############# utils ##########################

import Base: copy
function copy(state::ClusterWavefunction)
    tensors = copy([wf.itensor for wf in state.wavefunctions])
    indices = [wf.indices for wf in state.wavefunctions]
    return ClusterWavefunction([Wavefunction(t, inds) for (t, inds) in zip(tensors, indices)], state.indices)
end

function copy(gate::Gate)
    return Gate(gate.name, copy(gate.itensor))
end

import Base: +
function +(wf1::Wavefunction, wf2::Wavefunction)
    return Wavefunction(wf1.itensor + wf2.itensor, wf1.indices)
end


;
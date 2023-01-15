function get_line_topology(n_qubits)
    return [(ii, ii+1) for ii in 1:n_qubits-1]
end

function get_star_topology(n_qubits)
    return [(1, ii) for ii in 2:n_qubits]
end

function get_grid_topology(nx, ny)
    neighbors = []
    for y in 0:ny-1
        append!(neighbors, [(y*nx + ii, y*nx + ii+1) for ii in 1:nx-1])
    end
    for x in 1:nx
        append!(neighbors, [(x + ii*nx, x + (ii+1)*nx) for ii in 0:ny-2])
    end
    return neighbors
end

function get_rainbow_topology(n_qubits)
    return [(ii, n_qubits+1-ii) for ii in 1:floor(Int, n_qubits/2)]
end



function get_pauli_circuit(indices, n_layers::Int)

    n_qubits = length(indices)
    
    circuit::Vector{Gate} = []
    for nl in 1:n_layers
        for ii in 1:n_qubits
            append!(circuit, [get_x_gate(indices[ii])])
            append!(circuit, [get_z_gate(indices[ii])])
        end
        for ii in 1:n_qubits-1
            append!(circuit, [get_cnot_gate([indices[ii], indices[ii+1]])])
        end
    end
    return circuit
end


function get_ghz_circuit(indices, star_topology=false)
    n_qubits = length(indices)
    
    circuit::Vector{Gate} = []
    append!(circuit, [get_h_gate(indices[1])])
    if star_topology
        for ii in 2:n_qubits
            append!(circuit, [get_cnot_gate([indices[1], indices[ii]])])
        end
    else
        for ii in 1:n_qubits-1
            append!(circuit, [get_cnot_gate([indices[ii], indices[ii+1]])])
        end
    end
    
    return circuit
end

function get_random_circuit(indices, n_layers::Int, topology)
    
    n_qubits = length(indices)
    
    single_qubit_gate_choices = [get_sx_gate, get_sz_gate, get_t_gate]
    
    circuit::Vector{Gate} = []
    
    for ii in 1:n_qubits
        append!(circuit, [get_h_gate(indices[ii])])
    end
    for nl in 1:n_layers
        for ii in 1:n_qubits
            append!(circuit, [single_qubit_gate_choices[rand(1:3)](indices[ii])])
            append!(circuit, [single_qubit_gate_choices[rand(1:3)](indices[ii])])
        end
        for (ii, jj) in topology
            append!(circuit, [get_cnot_gate([indices[ii], indices[jj]])])
        end
    end
    return circuit
end

function get_random_circuit(indices, n_layers::Int)
    n_qubits = length(indices)
    return get_random_circuit(indices, n_layers, get_line_topology(n_qubits))
end


;
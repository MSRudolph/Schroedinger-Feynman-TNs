using ITensors

function get_x_gate(index::Index)
    indices = [index, index']
    matrix = [[0 1]; [1 0]]
    
    return Gate("X", ITensor(matrix, indices...))
end

function get_sx_gate(index::Index)
    indices = [index, index']
    matrix = [[1+1im 1-1im]; [1-1im 1+1im]] ./2
    
    return Gate("SX", ITensor(matrix, indices...))
end

function get_y_gate(index::Index)
    indices = [index, index']
    matrix = [[0 1im]; [-1im 0]]
    
    return Gate("Y", ITensor(matrix, indices...))
end

function get_z_gate(index::Index)
    indices = [index, index']
    matrix = [[1 0]; [0 -1]]
    
    return Gate("Z", ITensor(matrix, indices...))
end

function get_sz_gate(index::Index)
    indices = [index, index']
    matrix = [[1 0]; [0 1im]]
    
    return Gate("SZ", ITensor(matrix, indices...))
end

function get_h_gate(index::Index)
    indices = [index, index']
    matrix = 1/sqrt(2) .* [[1 1]; [1 -1]]
    
    return Gate("H", ITensor(matrix, indices...))
end

function get_t_gate(index::Index)
    indices = [index, index']
    matrix = [[1 0]; [0 exp(-π * 1im/ 4)]]
    
    return Gate("T", ITensor(matrix, indices...))
end


function get_cnot_gate(indices)
    indices = [indices..., indices'...]
    matrix = permutedims(reshape([
        [1 0 0 0];
        [0 1 0 0];
        [0 0 0 1];
        [0 0 1 0]
        ], (2, 2, 2, 2)), (4, 3, 2, 1))
    return Gate("CNOT", ITensor(matrix, indices...))
end

function get_custom_gate(matrix, indices)
    indices = [indices..., indices'...]
    return Gate("Gate", ITensor(matrix, indices...))
end

function get_custom_gate(matrix, index::Index)
    indices = [index, index']
    return Gate("Gate", ITensor(matrix, indices...))
end


function decompose(gate::Gate)
    if gate.name == "CNOT"
        indices = gate.itensor.tensor.inds
        p0_gate = get_custom_gate([[1 0]; [0 0]], indices[1])
        p1_gate = get_custom_gate([[0 0]; [0 1]], indices[1])
        x_gate = get_x_gate(indices[2])
        return [[p0_gate], [p1_gate; x_gate]]
    else
        return [gate]
    end
end
;

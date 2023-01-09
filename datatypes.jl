using ITensors


struct Gate 
    name::String
    itensor::ITensor
end

struct Wavefunction
    itensor::ITensor
    indices::Vector{Index}
end


struct ClusterWavefunction
    wavefunctions::Vector{Wavefunction}
    indices::Vector{Index}
end

;
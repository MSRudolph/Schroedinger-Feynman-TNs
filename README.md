# Schroedinger-Feynman-TNs
A repository with Julia implementations of several quantum circuit simulation methods: Typical statevector simulation, MPS simulation, and Feynman path-type simulation. The code is heavily based on the [Itensors.jl](https://github.com/ITensor/ITensors.jl) library.

For the slides used in the project presentation, see [here](https://github.com/MSRudolph/Schroedinger-Feynman-TNs/blob/main/schroedinger-feynman.pptx). Introductory figures as well as all benchmarks can be seen there.

Functionality and usability upgrades will follow.


## The Feynman Method
Because of the linearity of quantum mechanics, the fact that multi-qubit gates can be expressed as the sum of single-qubit gates operators implies that such entangling gates can be realized by applying the single-qubit operatos to different "replicas" of a wavefunction. Each wavefunction is then no longer normalized, but their sum represents the exact outcome of the entangling gate application without having ever applied a multi-qubit gate. \
The following example depicts the application of a CNOT gate as the sum of two single-qubit terms. In general, $n$-qubit gates can be decomposed into up to $2^n$ terms. This number is also called the _rank_ of the entangling gate.

<img src="https://github.com/MSRudolph/Schroedinger-Feynman-TNs/blob/main/figures/Feynman_technique.PNG" width="800">

A simulator based on the *Feynman* method scales exponentially in every multi-qubit gate that one applies. The *Schrödinger-Feynman* method finds a compromise of efficient statevector simulation and the number of replicas or paths. Here, the full wavefunction is split up into tensors of smaller wavefunctions. Any gates applied to qubits within one wavefunction is applied using statevector simulation. Gates applied to qubits across different wavefuntions are implemented via the Feynman method - in turn causing a doubling of the number of wavefunctions that are being simulated. While this method still scales exponentially in circuit depth, it now only scales exponentially in the number of gates crossing wavefunction borders. 

<img src="https://github.com/MSRudolph/Schroedinger-Feynman-TNs/blob/main/figures/S-F_technique.PNG" width="800">

This technique allows to simulate shallow quantum circuits acting on surprisingly large numbers of qubits.

## Example Benchmarks

### 1D Random Circuits

Results of exact statevector simulation and the Schrödinger-Feynman technique. The red and blue curves demonstrate that the choice of how to partition the wavefunction is very important for performance and feasibility. There is a trade-off between the exponential scaling in exactly simulating large wavefunctions and the number of cross-wavefunction gates.

<img src="https://github.com/MSRudolph/Schroedinger-Feynman-TNs/blob/main/figures/cluster_qubits_with_optimal_depth3_and_exact.png" width="700">

When directly comparing the Schrödinger-Feynmann method with optimal partitioning, we see that for lower number of qubits, it can be around one order of magnitude faster than MPS simulation.

<img src="https://github.com/MSRudolph/Schroedinger-Feynman-TNs/blob/main/figures/1d_mps_vs_cluster.png" width="700">

A study of circuit depth showcases that both the Schrödinger-Feynman method and *Matrix Product State* (MPS) simulation of quantum circuits scales exponentially, whereas statevector simulation scales only linearly in depth. If the full wavefunction however can be simualted, as in this 20 qubit example, the Schrödinger-Feynman method with optimized wavefunction size is always at least as fast as the statevector method. 

<img src="https://github.com/MSRudolph/Schroedinger-Feynman-TNs/blob/main/figures/scaling_with_depth.png" width="700">

### 2D Random Circuits

2D wavefunctions are significantly harder to simulate because there are more cross-wavefunction gates per layer. Thus, even though MPS simulation also scales exponentially with depth, Schrödinger-Feynman method carries a lot of redundant information and scales thus even worse. 

<img src="https://github.com/MSRudolph/Schroedinger-Feynman-TNs/blob/main/figures/cluster_and_mps_2d_exact.png" width="700">

### Fidelity = 0.5%

Finite fidelity can be achieved by dropping replicas (or paths) in the Feynman simulation. Since the many wavefunctions have exponentially similar weights, there is no clear way of choosing which to drop. For the same reason, one achieves a directly linear correspondence between the fidelity and the number of replicas. By dropping 99.5% of replicas, one achieves a fidelity of 0.5%. This is opposed to the MPS with approximates simulation not on the level of the gate (which is always the same), but on the level of the state. Thus, one can usually expect toachieve  drastic improvements in fidelity using the same memory resources as with the Schrödinger-Feynman method.

<img src="https://github.com/MSRudolph/Schroedinger-Feynman-TNs/blob/main/figures/truncation_efficiency.png" width="700">

Using the Schrödinger-Feynman method, one layer of the 1D circuit on 210 qubits can be simulated in one second.

<img src="https://github.com/MSRudolph/Schroedinger-Feynman-TNs/blob/main/figures/cluster_and_mps_1d_0005.png" width="700">

The key for being able to simulate 49 qubit 2D circuits more effectively is to be able to simulate 25 qubit using statevector simulation. However, this is likely not possible on your machine as the statevector simulation is rather naively implemented.

<img src="https://github.com/MSRudolph/Schroedinger-Feynman-TNs/blob/main/figures/cluster_and_mps_2d_0.005.png" width="700">

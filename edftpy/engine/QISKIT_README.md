# Qiskit Engine for eDFTpy

## Overview

The Qiskit engine (`EngineQiskit`) enables the use of Qiskit's Variational Quantum Eigensolver (VQE) to solve subsystems in embedded DFT calculations. This allows you to leverage quantum computing algorithms for quantum chemistry calculations within the eDFTpy framework.

## Requirements

To use the Qiskit engine, you need to install the following packages:

```bash
pip install qiskit
pip install qiskit-nature
pip install qiskit-algorithms
pip install qiskit-aer  # For statevector simulator
```

Optional (for real quantum hardware):
```bash
pip install qiskit-ibm-provider
```

## Configuration

### Basic Configuration

In your eDFTpy configuration file, set the calculator to `qiskit`:

```ini
[SUB0]
calculator = qiskit
technique = KS
cell.index = [0, 1, 2]
prefix = sub0_qiskit
```

### Qiskit-Specific Options

You can configure Qiskit-specific parameters in the `[SUB0]` section:

```ini
[SUB0]
# ... other options ...

# Qiskit options
qiskit.mapper = jordan_wigner    # Qubit mapping: jordan_wigner, parity, bravyi_kitaev
qiskit.ansatz = uccsd            # Ansatz type: uccsd, real_amplitudes
qiskit.optimizer = slsqp         # Classical optimizer: slsqp, spsa
qiskit.basis = sto-3g            # Basis set for PySCF (e.g., sto-3g, 6-31g, cc-pvdz)
qiskit.method = rhf               # Electronic structure method: rhf, uhf
# qiskit.backend = None           # Quantum backend (None = statevector simulator)
```

### Option Descriptions

- **mapper**: The qubit mapping strategy
  - `jordan_wigner`: Standard Jordan-Wigner transformation (default)
  - `parity`: Parity mapping (reduces qubit count for some systems)
  - `bravyi_kitaev`: Bravyi-Kitaev transformation

- **ansatz**: The quantum circuit ansatz
  - `uccsd`: Unitary Coupled Cluster Singles and Doubles (default, more accurate)
  - `real_amplitudes`: RealAmplitudes ansatz (simpler, fewer parameters)

- **optimizer**: Classical optimizer for VQE
  - `slsqp`: Sequential Least Squares Programming (default, gradient-based)
  - `spsa`: Simultaneous Perturbation Stochastic Approximation (noise-resistant)

- **basis**: Basis set for PySCF (standard quantum chemistry basis sets)
  - Examples: `sto-3g`, `6-31g`, `cc-pvdz`, `cc-pvtz`

- **method**: Electronic structure method
  - `rhf`: Restricted Hartree-Fock (closed-shell systems)
  - `uhf`: Unrestricted Hartree-Fock (open-shell systems)

- **backend**: Quantum backend
  - `None`: Uses statevector simulator (exact, classical simulation)
  - Or provide a Qiskit backend object for real quantum hardware

## Usage Example

See `examples/qiskit_example.ini` for a complete example configuration.

## Orbital Treatment and Wannier Functions

**Short answer: No, Wannierization is NOT required for the current Qiskit engine implementation.**

### How the Current Implementation Works

The current Qiskit engine implementation:

1. **Creates a standalone molecular calculation**: Each subsystem is treated as an isolated molecule
2. **Uses PySCF directly**: The engine creates a fresh PySCF calculation from the subsystem's atomic positions
3. **No pre-computed orbitals needed**: The engine doesn't use orbitals from the global system
4. **Molecular orbitals are generated**: PySCF generates molecular orbitals for the subsystem independently

### When Wannier Functions Might Be Useful

Wannier functions could be beneficial in these scenarios:

1. **Periodic systems**: If your global system is periodic and you want to use orbitals from the global calculation
2. **Better initial guess**: Using Wannier-localized orbitals from the global system as a starting point
3. **Orbital constraints**: If you want to constrain certain orbitals in the subsystem calculation
4. **Advanced embedding**: For more sophisticated embedding schemes that use orbital-based embedding

### Current Approach vs. Wannier-Based Approach

**Current Approach (Implemented)**:
```
Subsystem atoms → PySCF → Molecular orbitals → Qiskit VQE
```
- Simple and self-contained
- Works for molecular subsystems
- No dependency on global system orbitals

**Wannier-Based Approach (Future Enhancement)**:
```
Global system → Wannier functions → Localized orbitals → Subsystem → Qiskit VQE
```
- More integrated with global system
- Better for periodic systems
- Requires Wannier transformation step

### For Most Use Cases

For typical eDFTpy calculations with molecular subsystems, **you do NOT need Wannier functions**. The current implementation works well because:

- Subsystems are typically small molecular fragments
- Each subsystem is solved independently
- The embedding potential (not orbitals) connects subsystems
- PySCF generates appropriate molecular orbitals automatically

## Current Limitations

1. **Density Extraction**: The current implementation has a placeholder for density extraction from the quantum state. A full implementation would require:
   - Computing the 1-RDM (one-electron reduced density matrix) from the quantum state
   - Transforming to real-space density using basis functions
   - Mapping to the eDFTpy grid

2. **Forces and Stress**: Not yet implemented. Would require computing energy gradients.

3. **Large Systems**: VQE is currently most practical for small molecules (few atoms) due to:
   - Exponential scaling of quantum circuits
   - Limited qubit counts on current quantum hardware
   - Classical simulation overhead

4. **Embedding Potential**: The external embedding potential is stored but not yet fully integrated into the Hamiltonian construction. Currently, the subsystem is treated as isolated, and the embedding potential would need to be added to the Hamiltonian matrix elements.

5. **Orbital Integration**: The engine doesn't currently use orbitals from the global system. For advanced embedding schemes, you might want to:
   - Use Wannier functions to localize global orbitals
   - Use these as initial guess or constraints
   - This would require additional implementation

## How It Works

1. **Initialization**: The engine creates a PySCF driver from the molecular geometry
2. **Problem Setup**: Converts the electronic structure problem to a qubit Hamiltonian
3. **VQE Solution**: Uses VQE to find the ground state
4. **Energy Extraction**: Extracts the total energy from the VQE result

## Integration with eDFTpy

The Qiskit engine integrates seamlessly with eDFTpy's subsystem DFT framework:

- Can be used alongside other engines (QE, CASTEP, DFTpy)
- Supports embedding potential from other subsystems
- Follows the same interface as other engines

## Future Improvements

- Full density extraction from quantum states
- Force calculations via energy gradients
- Support for excited states (TDDFT)
- Better integration of embedding potentials
- Optimization for larger systems

## References

- Qiskit Nature: https://qiskit.org/ecosystem/nature/
- VQE Algorithm: https://qiskit.org/documentation/nature/stubs/qiskit_nature.second_q.algorithms.VQEUCCFactory.html
- eDFTpy Documentation: See main README


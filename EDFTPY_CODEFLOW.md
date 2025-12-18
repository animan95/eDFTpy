# eDFTpy Code Flow and Architecture

## Overview
eDFTpy (Embedded Density Functional Theory in Python) is a package for performing subsystem DFT calculations, where a large system is decomposed into smaller subsystems that are solved separately and then combined. It supports multiple quantum chemistry engines (Quantum ESPRESSO, CASTEP, DFTpy) and can perform various types of calculations (SCF, TDDFT, QM/MM).

---

## 1. Entry Point and Initialization

### 1.1 Command Line Interface
**File**: `edftpy/__main__.py` → `edftpy/cui/main.py`

```
python -m edftpy config.ini
```

**Flow**:
1. `main()` in `cui/main.py` parses command-line arguments
2. Routes to `run.main()` for calculations
3. Reads configuration file using `read_conf()`

### 1.2 Configuration Parsing
**File**: `edftpy/interface.py` → `edftpy/api/parse_config.py`

**Key Functions**:
- `read_conf()`: Reads INI/JSON config file
- `config2optimizer()`: Main initialization function
- `config2gsystem()`: Creates global system
- `config2driver()`: Creates subsystem drivers

**Process**:
```
config.ini → read_conf() → config2optimizer()
                                    ↓
                    ┌───────────────┴───────────────┐
                    ↓                               ↓
            config2gsystem()              config2driver() (for each SUB)
                    ↓                               ↓
            GlobalCell created          DriverKS/DriverOF created
                    ↓                               ↓
            TotalEvaluator              EngineQE/EngineDFTpy
            (functionals)               (wraps QE/DFTpy)
```

---

## 2. System Architecture

### 2.1 Global System (GSYSTEM)
**Class**: `GlobalCell` in `edftpy/subsystem/subcell.py`

**Components**:
- **Ions**: Atomic positions and cell
- **Grid**: Real-space grid for density representation
- **Density**: Total electron density (sum of all subsystems)
- **TotalEvaluator**: Functional evaluator for global system
  - Contains: HARTREE, XC, KE, PSEUDO functionals

**Purpose**: 
- Stores the complete system
- Computes embedding potentials
- Combines subsystem densities

### 2.2 Subsystems (SUB)
**Class**: `SubCell` in `edftpy/subsystem/subcell.py`

**Components**:
- **Ions**: Subset of atoms for this subsystem
- **Grid**: Subsystem grid (may differ from global)
- **Density**: Subsystem electron density
- **EmbedEvaluator**: Functional evaluator for subsystem
  - Computes embedding potential from global system

**Purpose**:
- Represents a fragment of the total system
- Solved independently with embedding potential

### 2.3 Drivers
**Classes**: `DriverKS`, `DriverOF`, `DriverMM` in `edftpy/engine/driver.py`

**DriverKS** (Kohn-Sham):
- Uses external engines (QE, CASTEP) via `EngineQE`, `EngineCastep`
- Performs SCF/TDDFT calculations
- Manages density mixing

**DriverOF** (Orbital-Free):
- Uses DFTpy directly
- Optimizes density without orbitals

**DriverMM** (Molecular Mechanics):
- Classical force field calculations

---

## 3. Main Calculation Flow

### 3.1 Optimization Loop
**File**: `edftpy/optimizer.py` - `Optimization.optimize()`

```
┌─────────────────────────────────────────┐
│ 1. Initialization                       │
│    - Build global system                │
│    - Create subsystems                  │
│    - Initialize drivers                 │
└─────────────────────────────────────────┘
                    ↓
┌─────────────────────────────────────────┐
│ 2. Initial Density Setup                │
│    - Sum subsystem densities            │
│    - Initialize global density          │
└─────────────────────────────────────────┘
                    ↓
┌─────────────────────────────────────────┐
│ 3. Optimization Loop (maxiter times)    │
│                                         │
│    ┌───────────────────────────────┐   │
│    │ 3a. Set Global Potential      │   │
│    │     - Compute embedding pot   │   │
│    │     - Distribute to subs      │   │
│    └───────────────────────────────┘   │
│                  ↓                      │
│    ┌───────────────────────────────┐   │
│    │ 3b. Update Each Subsystem      │   │
│    │     - Set external potential   │   │
│    │     - Run SCF/TDDFT           │   │
│    │     - Get new density         │   │
│    └───────────────────────────────┘   │
│                  ↓                      │
│    ┌───────────────────────────────┐   │
│    │ 3c. Update Global Density     │   │
│    │     - Sum subsystem densities  │   │
│    │     - Mix density              │   │
│    └───────────────────────────────┘   │
│                  ↓                      │
│    ┌───────────────────────────────┐   │
│    │ 3d. Check Convergence         │   │
│    │     - Energy change           │   │
│    │     - Density change         │   │
│    └───────────────────────────────┘   │
└─────────────────────────────────────────┘
                    ↓
┌─────────────────────────────────────────┐
│ 4. Finalize                             │
│    - Calculate forces/stress            │
│    - Write output files                 │
└─────────────────────────────────────────┘
```

### 3.2 Step-by-Step Details

#### Step 3a: Set Global Potential
**Function**: `Optimization.set_global_potential_sdft()`

**Process**:
1. Compute embedding potential from global density:
   ```python
   embed_potential = TotalEvaluator.get_embed_potential(global_density)
   ```
   - This includes: Hartree, XC, KE, Pseudo potentials
   - Excludes subsystem-specific terms

2. Distribute to each subsystem:
   ```python
   gsystem.sub_value(embed_potential, driver.evaluator.global_potential, isub=i)
   ```
   - Extracts subsystem region from global potential
   - Stores in `driver.evaluator.global_potential`

**Mathematical Formulation**:
```
v_embed^I = v_H[ρ_total] + v_XC[ρ_total] + v_KE[ρ_total] 
          - Σ_J≠I (v_H[ρ_J] + v_XC[ρ_J] + v_KE[ρ_J])
```

#### Step 3b: Update Subsystem
**Function**: `DriverKS.__call__()` → `Driver.compute()`

**For KS Drivers (QE, CASTEP)**:
1. Set external potential:
   ```python
   engine.set_extpot(embedding_potential)
   ```

2. Run SCF:
   ```python
   engine.scf()  # Diagonalizes Hamiltonian
   ```

3. Get density:
   ```python
   engine.get_rho(density)
   ```

**For OF Drivers (DFTpy)**:
1. Optimize density directly:
   ```python
   evaluator_of.compute(density, calcType=['E','V'])
   # Minimizes energy functional
   ```

#### Step 3c: Update Global Density
**Function**: `Optimization.update_density()`

**Process**:
1. Gather subsystem densities:
   ```python
   for driver in drivers:
       gsystem.update_density(driver.density, isub=i)
   ```

2. Mix density (if mixer enabled):
   ```python
   density_new = mixer(density_prev, density_new)
   ```

3. Update global density:
   ```python
   gsystem.density = sum of all subsystem densities
   ```

#### Step 3d: Check Convergence
**Functions**: 
- `check_converge_energy()`: Energy change < econv
- `check_converge_potential()`: Density change < pconv

---

## 4. Engine Interface (Quantum ESPRESSO Example)

### 4.1 EngineQE Class
**File**: `edftpy/engine/engine_qe.py`

**Purpose**: Wraps QEpy (Python interface to Quantum ESPRESSO)

**Key Methods**:
- `initial()`: Initialize QE calculation
  - Creates input file
  - Starts QE driver
  - Sets up embedding

- `scf()`: Run SCF calculation
  - Calls `driver.diagonalize()`
  - Solves Kohn-Sham equations

- `set_extpot()`: Set embedding potential
  - Passes external potential to QE

- `get_rho()`: Get electron density
  - Extracts density from QE

- `get_energy()`: Get total energy
  - Returns QE total energy

### 4.2 Driver-Engine Interaction
```
DriverKS
    ↓
    engine.initial()          # Setup QE
    ↓
    engine.set_extpot()      # Set embedding
    ↓
    engine.scf()             # Run SCF
    ↓
    engine.get_rho()         # Get density
    ↓
    engine.get_energy()      # Get energy
```

---

## 5. Evaluator System

### 5.1 Evaluator Hierarchy
**File**: `edftpy/evaluator.py`

**Classes**:
1. **Evaluator**: Base class
   - Computes functionals (HARTREE, XC, KE, PSEUDO)
   - Returns energy and potential

2. **EmbedEvaluator**: For subsystems
   - Computes embedding potential
   - Handles global potential

3. **TotalEvaluator**: For global system
   - Computes total functional
   - Generates embedding potential

### 5.2 Functional Computation
```python
# Example: Compute Hartree energy and potential
hartree = Hartree()
result = hartree(density, calcType=['E','V'])
# result.energy = Hartree energy
# result.potential = Hartree potential
```

---

## 6. Density Mixing

### 6.1 Mixer
**Class**: `Mixer` in `edftpy/mixer/`

**Purpose**: Stabilize SCF convergence

**Types**:
- **Pulay**: Uses history of previous densities
- **Linear**: Simple linear mixing

**Usage**:
```python
density_new = mixer(density_prev, density_new, coef=0.7)
```

---

## 7. Parallelization

### 7.1 MPI Structure
**File**: `edftpy/mpi/mpi.py`

**GraphTopo**:
- Manages MPI communicators
- Distributes processors to subsystems
- Handles communication between subsystems

**Structure**:
```
COMM_WORLD (all processors)
    ↓
    ├─ comm (global system)
    │
    └─ comm_sub (subsystems)
        ├─ SUB0: processors 0-3
        ├─ SUB1: processors 4-7
        └─ SUB2: processors 8-11
```

---

## 8. Output and Post-Processing

### 8.1 Output Files
**Function**: `conf2output()` in `edftpy/interface.py`

**Outputs**:
- Density files (`.xsf`, `.cube`)
- Embedding potential
- Forces
- Stress tensor
- Energy breakdown

### 8.2 Properties
**File**: `edftpy/properties/`

- **Forces**: `get_total_forces()`
- **Stress**: `get_total_stress()`
- **Electrostatic potential**: `get_electrostatic_potential()`

---

## 9. Special Calculation Types

### 9.1 TDDFT (Time-Dependent DFT)
**File**: `edftpy/tddft/optical.py`

**Class**: `MoleculeOpticalAbsorption`

**Flow**:
1. Run SCF to get ground state
2. Switch to TDDFT mode
3. Compute optical absorption

### 9.2 QM/MM
**Mode**: `sdft='qmmm'`

**Components**:
- QM subsystem: Quantum mechanics (DFT)
- MM subsystem: Molecular mechanics (classical)

**Interaction**: 
- QM density affects MM
- MM potential affects QM

### 9.3 PDFT (Partition DFT)
**Mode**: `sdft='pdft'`

**Difference from SDFT**:
- Uses different embedding potential formulation
- Better for charge transfer

---

## 10. Key Data Structures

### 10.1 Density
**Type**: `Field` (from DFTpy)

**Properties**:
- `grid`: Real-space grid
- `rank`: Spin dimension (1 or 2)
- `integral()`: Total charge

### 10.2 Grid
**Type**: `Grid` (from DFTpy)

**Properties**:
- `nr`: Grid points per dimension
- `lattice`: Unit cell vectors
- `dV`: Volume element

### 10.3 Ions
**Type**: `Ions` (from DFTpy)

**Properties**:
- `positions`: Atomic positions
- `cell`: Unit cell
- `symbols`: Element symbols

---

## 11. Configuration File Structure

```ini
[GSYSTEM]
# Global system settings
grid = [64, 64, 64]
decompose = distance
rcut = 3.0

[SUB0]
# Subsystem 0
calculator = qe
technique = KS
cell.index = [0, 1, 2, 3]  # Atom indices
prefix = sub0

[SUB1]
# Subsystem 1
calculator = dftpy
technique = OF
cell.index = [4, 5, 6]
prefix = sub1

[PP]
# Pseudopotentials
H = h_pbe.UPF
O = o_pbe.UPF

[JOB]
calctype = [Energy, Force]
```

---

## 12. Summary of Key Concepts

1. **Subsystem Decomposition**: Large system → smaller subsystems
2. **Embedding Potential**: Potential from other subsystems
3. **Iterative Solution**: 
   - Solve subsystems with embedding
   - Combine densities
   - Update embedding
   - Repeat until convergence
4. **Engine Abstraction**: Unified interface to different QM codes
5. **Parallel Execution**: Each subsystem can run on separate processors

---

## 13. Code Flow Summary

```
User Input (config.ini)
    ↓
Parse Config
    ↓
Create Global System (GSYSTEM)
    ↓
Decompose into Subsystems (SUB0, SUB1, ...)
    ↓
Create Drivers (DriverKS, DriverOF)
    ↓
Initialize Engines (EngineQE, EngineDFTpy)
    ↓
┌─────────────────────────────────────┐
│ OPTIMIZATION LOOP                   │
│                                     │
│ 1. Compute Embedding Potential     │
│    (from global density)            │
│                                     │
│ 2. For each subsystem:             │
│    a. Set embedding potential      │
│    b. Run SCF/TDDFT                │
│    c. Get new density              │
│                                     │
│ 3. Combine subsystem densities     │
│    → Update global density          │
│                                     │
│ 4. Mix density (if needed)         │
│                                     │
│ 5. Check convergence               │
│    - Energy change                 │
│    - Density change               │
│                                     │
│ Repeat until converged             │
└─────────────────────────────────────┘
    ↓
Calculate Properties (Forces, Stress)
    ↓
Write Output Files
    ↓
Done
```

---

This architecture allows eDFTpy to:
- Handle large systems by decomposition
- Use different QM engines for different subsystems
- Parallelize across subsystems
- Perform various types of calculations (SCF, TDDFT, QM/MM)


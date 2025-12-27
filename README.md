# eDFTpy: Embedding, Multiscale, and Hybrid Quantum–Classical Modeling in Python  
A creation of PRG

`eDFTpy` is a Python framework for **embedding-based and multiscale electronic-structure simulations**, combining **Kohn–Sham Density Functional Theory (KS-DFT)** and **orbital-free DFT (OF-DFT)** within a unified workflow. Weak intermolecular interactions are treated at the OF-DFT level, while strong intramolecular interactions are handled using KS-DFT, enabling efficient and scalable simulations of complex systems.

Building on this classical foundation, `eDFTpy` also provides a **hybrid quantum–classical execution layer** that enables **variational quantum algorithms (e.g., VQE, EOM-VQE)** to be integrated into embedding and first-principles workflows. This allows selected subsystems to be mapped to qubit Hamiltonians and executed on **quantum simulators or real quantum hardware**, while the remainder of the system is treated classically.

The code is built on top of [`DFTpy`](http://dftpy.rutgers.edu) and [`QEpy`](https://gitlab.com/shaoxc/qepy), and is developed by the [PRG](https://sites.rutgers.edu/prg/) at [Rutgers University–Newark](http://sasn.rutgers.edu).

---

## Key Capabilities

- Embedding-based multiscale simulations using **KS-DFT and OF-DFT**
- Seamless integration with **Quantum ESPRESSO** via QEpy
- **Hybrid quantum–classical workflows**, enabling:
  - Qubit Hamiltonian generation from first-principles calculations  
  - Execution of VQE and related algorithms  
  - Support for local simulators and quantum hardware backends
- Modular, extensible architecture designed for **scientific software development and experimentation**

---

## Requirements

- Python 3.6 or later  
- [DFTpy](https://gitlab.com/pavanello-research-group/dftpy) (latest)  
- [NumPy](https://numpy.org/doc/stable)  
- [SciPy](https://docs.scipy.org/doc/scipy/reference)  
- [ASE](http://wiki.fysik.dtu.dk/ase)

### Optional (highly recommended)

- [Libxc](https://gitlab.com/libxc/libxc)  
- [f90wrap](https://github.com/jameskermode/f90wrap) (latest)

### Optional (for hybrid quantum–classical workflows)

- [Qiskit](https://qiskit.org) (for variational quantum algorithm execution)

---

## Documentation

- [`eDFTpy` manual](http://edftpy.rutgers.edu)

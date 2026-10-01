# Monte Carlo Simulation of a Lennard–Jones Fluid

Computational physics project developed during my M.Sc. in Physics at **Sapienza University of Rome**.

This project implements Monte Carlo simulations of a Lennard–Jones fluid in different thermodynamic ensembles and compares independent approaches for estimating its thermodynamic properties and equation of state.

The main methods implemented are:

- **Isothermal-isobaric Monte Carlo (NPT)**
- **Grand-canonical Monte Carlo (μVT)**
- **Successive Umbrella Sampling (SUS)**
- **Chemical-potential reweighting**
- Comparison of equations of state obtained using different simulation approaches

---

## Physical model

The simulated system consists of particles interacting through the Lennard–Jones potential:

```math
U(r) =
4\epsilon
\left[
\left(\frac{\sigma}{r}\right)^{12}
-
\left(\frac{\sigma}{r}\right)^6
\right].
```

The simulations are performed using reduced Lennard–Jones units.

A cubic simulation box with **periodic boundary conditions** is used, together with the **minimum-image convention** for particle interactions.

---

## NPT Monte Carlo

The isothermal-isobaric implementation samples configurations at constant:

- particle number \(N\)
- pressure \(P\)
- temperature \(T\)

The simulation includes:

- single-particle displacement moves
- isotropic volume-change moves
- Metropolis acceptance criteria
- periodic boundary conditions
- calculation of thermodynamic observables
- estimation of density as a function of pressure
- equation-of-state analysis

The NPT simulations provide an independent estimate of the thermodynamic behaviour of the Lennard–Jones fluid.

---

## Grand-canonical Monte Carlo

The grand-canonical implementation samples configurations at fixed:

- chemical potential μ
- volume \(V\)
- temperature \(T\)

The number of particles is allowed to fluctuate during the simulation.

Implemented Monte Carlo moves include:

- particle displacement
- particle insertion
- particle deletion

The simulations generate particle-number distributions \(P(N)\), which can be used to study density fluctuations and thermodynamic behaviour in the grand-canonical ensemble.

---

## Successive Umbrella Sampling

**Successive Umbrella Sampling (SUS)** is used to improve the sampling of the particle-number distribution, particularly in regions that are difficult to explore using conventional grand-canonical Monte Carlo.

The particle-number range is divided into overlapping sampling windows.

The relative probabilities obtained in the different windows are then combined to reconstruct the global particle-number distribution.

This approach makes it possible to obtain a more detailed estimate of the statistical weight associated with different particle numbers.

---

## Chemical-potential reweighting

Once the particle-number probability distribution has been reconstructed, it can be reweighted to nearby values of the chemical potential.

This allows thermodynamic observables to be estimated over a range of chemical potentials without performing a completely independent simulation at every point.

The reconstructed distributions can therefore be used to study quantities such as:

- average particle number
- density
- particle-number fluctuations
- equation-of-state behaviour

---

## Comparison of simulation methods

A central part of the project is the comparison between thermodynamic results obtained from:

1. direct **NPT Monte Carlo simulations**
2. **μVT simulations combined with Successive Umbrella Sampling and reweighting**

The comparison provides an internal consistency check between two independent statistical-mechanics approaches.

---

## Repository structure

```text
NPT/
├── source and header files for NPT simulations
├── analysis scripts
├── simulation outputs
└── figures and equation-of-state results

muVT/
├── source and header files for grand-canonical simulations
├── Successive Umbrella Sampling implementation
├── analysis and reweighting scripts
├── simulation outputs
└── particle-number distributions

confronto/
├── data used for comparison between NPT and μVT/SUS
└── plots comparing the resulting equations of state
```

Some directories also contain intermediate simulation outputs, parameter files and figures used during development and analysis.

---

## Implementation

The main simulation code is written in **C++**.

The implementation includes custom components for:

- particle representation
- simulation-box geometry
- periodic boundary conditions
- minimum-image distances
- random-number generation
- Monte Carlo moves
- thermodynamic parameters
- data collection and histogram generation

Additional data analysis and visualization were performed using **Python**.

---

## Example analyses

The project contains examples of:

- NPT equation-of-state calculations
- grand-canonical particle-number distributions
- Successive Umbrella Sampling
- reconstruction of statistical weights
- chemical-potential reweighting
- comparison of independently estimated equations of state
- thermodynamic observable analysis
- simulation time-series analysis

Representative plots and numerical results are included in the corresponding project directories.

---

## Scientific background

Monte Carlo methods provide a numerical approach for sampling equilibrium probability distributions in statistical mechanics.

Different ensembles describe systems under different thermodynamic constraints. Comparing results obtained in the NPT and grand-canonical ensembles provides a useful way to test numerical implementations and investigate the consistency of thermodynamic predictions.

Successive Umbrella Sampling extends conventional Monte Carlo sampling by enabling the exploration of particle-number states that may otherwise be poorly sampled.

---

## Academic context

This project was developed as part of the **Physics of Liquids** course during my M.Sc. in Physics at Sapienza University of Rome.

The work provided practical experience in:

- statistical mechanics
- Monte Carlo methods
- numerical simulation
- scientific programming in C++
- probability distributions and statistical sampling
- thermodynamic ensembles
- quantitative data analysis
- Python-based analysis and visualization
- validation and comparison of independent computational methods

---

## Technologies

- **C++**
- **Python**
- Monte Carlo simulation
- Numerical data analysis
- Scientific visualization

---

## Author

**Marco Adinolfi Falcone**  
M.Sc. in Physics — Sapienza University of Rome

Research interests: MR physics, quantitative imaging, neuroimaging and scientific computing.

GitHub: [MarcoAdinolfiFalcone](https://github.com/MarcoAdinolfiFalcone)

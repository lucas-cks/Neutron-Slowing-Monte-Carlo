# Monte Carlo Simulation of Neutron Slowing Down

**Author:** Ching Kai Sing, Lucas 

## 1. Overview
[cite_start]This project implements a high-fidelity Monte Carlo simulation of neutron moderation in an infinite, capturing medium[cite: 1]. [cite_start]It is based on the classic 1956 research by **Coveyou, Bate, and Osborn** at the Oak Ridge National Laboratory[cite: 2]. [cite_start]The simulation explores the effects of moderator temperature on neutron flux spectra, validating results against **Wigner-Wilkins (1944)** theory[cite: 2].

## 2. Physics Background
The simulation models neutron thermalization with the following parameters:
* [cite_start]**Target Motion:** Moderator atoms follow a Maxwell-Boltzmann distribution[cite: 4].
* [cite_start]**Cross-sections:** $1/v$ absorption cross-section is assumed[cite: 4].
* [cite_start]**Scattering:** Elastic scattering with temperature-dependent target velocity sampling[cite: 4, 9].
* **Reproduction of Key Phenomena:**
    * [cite_start]**Spectral Hardening:** The thermal peak shifts to higher energies as absorption ($K$) increases[cite: 7].
    * [cite_start]**1/v Tail:** High-speed flux follows the theoretical $1/v$ behavior[cite: 7].
    * [cite_start]**Temperature Shift:** Deviation from the ideal Maxwell-Boltzmann distribution[cite: 7].

## 3. Implementation Details
* [cite_start]**Language:** C (for high-performance simulation) and Python (for data analysis)[cite: 3, 6].
* [cite_start]**Scale:** 1,000,000 neutron histories per case to ensure low statistical noise[cite: 9].
* [cite_start]**Sampling:** Implements rejection sampling for target velocities as described in the 1956 paper[cite: 9].
* [cite_start]**Normalization:** Tallies are correctly normalized to flux for comparison with analytical models[cite: 9].

## 4. Repository Structure
* [cite_start]`monte_carlo_neutron_slowing.c`: Main simulation engine in C[cite: 3].
* [cite_start]`plot.py`: Python script for visualization and statistical analysis[cite: 3].
* [cite_start]`data/`: CSV files containing simulation results for Hydrogen, Carbon, and High Absorption cases[cite: 3].
* [cite_start]`results/`: Generated plots showing flux comparisons and hardening factors[cite: 3].
* [cite_start]`docs/`: Original reference papers by Coveyou et al. and Wigner-Wilkins[cite: 3].

## 5. Usage
### Compilation
```bash
gcc -o monte_carlo_neutron_slowing monte_carlo_neutron_slowing.c -lm -O3
```

### Execution
```bash
./monte_carlo
python plot.py
```

## 6. Key Results & Validation
[cite_start]The simulation reproduces the theoretical relationship between moderator temperature ($T_m$) and effective neutron temperature ($T_e$)[cite: 8]:
$$T_m/T_e = 1 + 1.11 \times A \times K$$

| Case | Moderator | A | K | Peak Position ($y$) |
| :--- | :--- | :--- | :--- | :--- |
| 1 | Hydrogen | 1 | 0.18 | $\approx 1.094$ |
| 2 | Carbon | 12 | 0.10 | $\approx 1.406$ |
| 3 | High Absorption | 1 | 0.36 | $\approx 1.219$ |


## 7. References
1. Coveyou, R. R., Bate, R. R., & Osborn, R. K. (1956). "Effect of Moderator Temperature upon Neutron Flux in Infinite, Capturing Medium". [cite_start]*Journal of Nuclear Energy*. [cite: 11]
2. Wigner, E. P., & Wilkins, J. E. (1944). "Effect of the Temperature of the Moderator on the Velocity Distribution of Neutrons". [cite_start]*AECD-2275*. [cite: 12]

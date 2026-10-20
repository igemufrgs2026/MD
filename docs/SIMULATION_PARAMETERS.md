# Simulation parameters

| Parameter | Setting |
| --- | --- |
| GROMACS | 2023.3-Ubuntu_2023.3_1ubuntu3 (study-reported version) |
| Force field | CHARMM36, November 2018 (study-reported release) |
| Integration timestep | 0.002 ps (2 fs) |
| Continuation block | 2,500,000 steps = 5 ns |
| Initial production segment | 1,650,000 steps = 3.3 ns, advancing from 1.7 to 5 ns |
| Study simulation endpoint | 500 ns on the simulation time axis |
| Equilibration | Eight stages totaling 1.7 ns |
| Production temperature | 310 K |
| Thermostat | V-rescale; protein and non-protein groups; tau_t = 0.1 ps for each |
| Production pressure | 1 bar |
| Barostat | Isotropic Parrinello–Rahman; tau_p = 2 ps |
| Compressibility | 4.5 × 10⁻⁵ bar⁻¹ |
| Solvent | TIP3P water |
| Ions | 0.15 M NaCl and charge neutralization |
| Box | Cubic; study-reported side length 10.22 nm |
| Minimization | Restrained steepest descent, unrestrained steepest descent, conjugate gradients |
| Constraints | All bonds during MD |
| Electrostatics | PME |
| Coulomb and van der Waals cutoffs | 1.2 nm in the production templates |
| Coordinate and velocity output | Every 5,000 steps = 10 ps |
| Compressed trajectory output | Every 2,500 steps = 5 ps |
| Energy and log output | Every 1,000 steps = 2 ps |

## Equilibration schedule

| Stage | Ensemble / temperature | Duration |
| --- | --- | --- |
| EQ1 | NVT, 100 K | 300 ps |
| EQ2 | NVT, 200 K | 100 ps |
| EQ3 | NVT, 280 K | 100 ps |
| EQ4 | NPT, 280 K, Berendsen | 200 ps |
| EQ5 | NPT, 280 K, Berendsen | 200 ps |
| EQ6 | NPT, 280 K, Berendsen | 200 ps |
| EQ7 | NPT, 290 K, Parrinello–Rahman | 200 ps |
| EQ8 | NPT, 310 K, Parrinello–Rahman | 400 ps |

## Interpreting the files

The `nsteps` setting in `070_md_full_5ns-step.mdp` is a block length, not the complete trajectory duration. The main protein–ligand pipeline sets `ENDTIME_NS=500`; the separate preserved preparation pipeline sets `ENDTIME_NS=10` and requires configuration for a longer run.

Use active MDP assignments rather than historical values in comments. Command-line overrides in the selected script can also affect the effective run settings. The original execution logs are not included here, so this collection documents the configured workflow rather than independently verifying every executed run.

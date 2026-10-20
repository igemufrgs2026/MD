# CadR molecular dynamics — iGEM UFRGS Brazil 2026

Input parameters and shell scripts for the molecular dynamics workflow discussed on our [Model wiki page](https://2026.igem.wiki/ufrgs-brazil/model).

## Documentation

- [Usage and required inputs](docs/USAGE.md)
- [Simulation parameters](docs/SIMULATION_PARAMETERS.md)
- [Citation information](CITATION.cff)
- Questions: [igemufrgs@gmail.com](mailto:igemufrgs@gmail.com)

## Files

- `MDPs_novos/005_em_pr.mdp`, `010_em_full_steep.mdp` and `015_em_full_cg.mdp`: energy minimization.
- `MDPs_novos/020_md_eq_1.mdp` through `055_md_eq_8.mdp`: eight equilibration stages.
- `MDPs_novos/060_md_full_005ns.mdp`: initial production segment, from 1.7 to 5 ns.
- `MDPs_novos/070_md_full_5ns-step.mdp`: continued production in 5 ns blocks.
- `MDPs_novos/run_md_corrigido.sh`: preparation and simulation pipeline; this copy has `ENDTIME_NS=10` and local paths that must be configured for a new run.
- `MD_script_CGenFF.sh`: protein–ligand preparation and simulation pipeline, configured with `ENDTIME_NS=500`.
- `concat_pbc_fit.sh`: trajectory processing.
- `analyze_md_fixed.sh`: protein–ligand trajectory analysis.

## Simulation settings

The study used GROMACS 2023.3-Ubuntu_2023.3_1ubuntu3 and the November 2018 CHARMM36 release, with TIP3P water. The reported simulation box was cubic, with a side length of 10.22 nm. Preparation scripts add 0.15 M NaCl and neutralize the system.

The integration timestep is 0.002 ps. The continued production template uses 2,500,000 steps **per 5 ns block**, not for the entire run. The reported trajectories extend to 500 ns on the simulation time axis, including 1.7 ns of equilibration. The first production segment has 1,650,000 steps and extends from 1.7 to 5 ns.

Production uses V-rescale temperature coupling at 310 K and isotropic Parrinello–Rahman pressure coupling at 1 bar. Coordinates and velocities are written every 5,000 steps (10 ps), compressed coordinates every 2,500 steps (5 ps), and energies/logs every 1,000 steps (2 ps).

Equilibration comprises NVT heating, NPT stabilization with Berendsen pressure coupling and progressive restraint release, followed by NPT stages with Parrinello–Rahman coupling.

## Using the scripts

Inspect each script before running it. Configure local paths, input structures, index groups, GROMACS executable and available hardware for your environment. The shell scripts require Bash; optional stages require the external tools referenced within each script. The production length must be checked in the selected pipeline.

The original study files are preserved without rewriting their settings. Some comments retain earlier parameter values: the active MDP assignments determine the configured values. These files are workflow inputs and scripts, not a complete archive of trajectories, checkpoints, structures or execution logs.

## Citation

Please cite the repository and the associated wiki page when referring to this workflow. GitHub can display citation metadata from `CITATION.cff`. Cite GROMACS, the force field and any other software separately when using those methods in your own work.

## Repository scope

This repository accompanies the UFRGS Brazil 2026 iGEM project. It contains study-specific scripts and input parameters. It does not contain measured binding affinities or experimental validation of the simulated interactions.

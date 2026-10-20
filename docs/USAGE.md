# Usage and required inputs

## Requirements

- Linux or another compatible Bash environment with GNU command-line tools.
- GROMACS; the study used `2023.3-Ubuntu_2023.3_1ubuntu3`. Scripts invoke `gmx_mpi` by default.
- Python 3 for embedded topology-processing and analysis helpers; these helpers use the Python standard library.
- CHARMM36 (November 2018 release), with the required water and ion topology files.
- Compatible GPU support if retaining the GPU flags in the simulation scripts.

GROMACS, force-field files, protein structures, ligand structures and ligand parameters must be installed or prepared separately. The repository does not bundle these dependencies.

## Choose a workflow

### Protein preparation and simulation

`MDPs_novos/run_md_corrigido.sh` prepares topology, solvent and ions, performs minimization and equilibration, and continues production in blocks. Set `MDP_DIR`, `CHARMM36_DIR`, input structures, group selections and `ENDTIME_NS` before execution. This preserved copy has a 10 ns endpoint; the study trajectories were extended to 500 ns.

### Protein–ligand preparation and simulation

`MD_script_CGenFF.sh` is configured with a 500 ns endpoint. It expects a prepared protein PDB, a ligand PDB, a compatible ligand topology, the MDP directory and the force-field directory. Its preserved input filenames are `protein_a8119c2bf1_prep_fixed.pdb`, `bestranking_gmx.pdb` and `bestranking_gmx.top`. Update these and the absolute local paths for the intended system. The script contains checks tied to the original system, including atom counts, which must also match the new input.

## Before launching

1. Use a separate working directory for each simulation or replicate.
2. Copy the selected script and configure its input paths and settings.
3. Inspect the active assignments in every MDP file. Some comments describe previous settings.
4. Confirm atom order, topology compatibility, index groups and installed software.
5. Check syntax with `bash -n script-name.sh` before executing the configured script with Bash.

The scripts prepare systems, create output files and may replace existing outputs. Run them in a dedicated simulation directory.

## Trajectory processing

`concat_pbc_fit.sh` concatenates matching trajectory segments, removes periodic-boundary jumps and fits rotation/translation using the protein by default. It expects `target_trajXTC_MD_*ns.xtc` files and a reference TPR. Settings can be changed through `GMX`, `REF`, `FIT_GROUP` and `NOWATER_GROUP` environment variables. Check whether files are separate segments or cumulative trajectories before concatenating them.

## Analysis

`analyze_md_fixed.sh` uses the final trajectory directly. Its default inputs include `targetmdrun_MD_500ns.tpr`, `target_trajXTC_MD_500ns.xtc`, `energy_MD_500ns.edr` and `target_BOX_MD_500ns.gro`. Override `GMX`, `REF_TPR`, `TRAJ`, `ENERGY`, `STRUCTURE` and `LIGAND_RESNAME` where necessary. Outputs are organized under `analysis/`.

The script checks input trajectories and records software metadata. Read the generated logs and inspect the selections before interpreting the results.

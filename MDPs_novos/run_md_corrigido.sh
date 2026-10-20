#!/usr/bin/env bash
set -euo pipefail

############################################################
# Pipeline MD corrigido - GROMACS + CHARMM36
#
# MODIFICACOES PRINCIPAIS NO SCRIPT:
#
# [A] Sequencia de equilibracao completamente refeita:
#     - EQ1-EQ3: NVT puro (sem barostato)
#       Aquecimento gradual: 100K -> 200K -> 280K
#     - EQ4-EQ6: NPT com Berendsen
#       Soltura gradual dos restraints: 100 -> 10 -> 0
#       Berendsen para estabilizar volume antes do PR
#     - EQ7-EQ8: NPT com Parrinello-Rahman
#       Aquecimento final: 290K -> 310K
#       PR so entra depois que volume ja esta equilibrado
#
# [B] update_posre chamado nos passos corretos:
#     - posre 5000 ja e o default dos arquivos .itp originais
#     - posre 1000 -> EQ3 (NVT, ainda aquecendo)
#     - posre 100  -> EQ4 (inicio NPT Berendsen)
#     - posre 10   -> EQ5
#     - posre 0    -> EQ6 (sem restraints, ainda Berendsen)
#
# [C] Flags GPU adicionadas tambem na equilibracao (EQ4+)
#     MOTIVO: os passos de equilibracao NPT sao mais longos
#     agora (200 ps cada). Vale usar GPU para acelerar.
#     NVT nao usa GPU pois sao passos curtos e sem barostato.
#
# [D] tinit do loop de producao corrigido:
#     O offset agora leva em conta o tempo total de
#     equilibracao (~1700 ps) para manter continuidade.
############################################################

########## USER SETTINGS ###################################
ENDTIME_NS=10
MDP_DIR="/home/administrator/Documents/MD_FILES_updated/mdps_corrigidos/"
CHARMM36_DIR="/home/administrator/Documents/MD_FILES_updated/md_files/charmm36.ff"
WATERMODEL="tip3p"
FFNAME="charmm36"
PME_THREADS="-nb gpu -bonded gpu -pme gpu -pin on -nstlist 400"
#-nb gpu	Non-bonded interactions on GPU
#-bonded gpu	Bonded interactions on GPU
#-pme gpu	Particle-Mesh Ewald on GPU
#-update gpu	Coordinate & velocity updates on GPU (fastest)
#-nstlist 400	Rebuilds neighbor list every 400 steps (reduces CPU load)
############################################################

INPUT_PDB=$(ls *.pdb 2>/dev/null | head -n1 || true)
INPUT_NDX=$(ls *.ndx 2>/dev/null | head -n1 || true)

if [[ -z "${INPUT_PDB}" ]]; then
    echo "ERROR: no .pdb file found in current directory."
    exit 1
fi

echo "Using PDB: ${INPUT_PDB}"
[[ -n "${INPUT_NDX}" ]] && echo "Using existing index: ${INPUT_NDX}"

STARTTIME=$(date "+%H:%M:%S-%D")
echo "Run started at ${STARTTIME}"

############################################################
# 2) topology + box + solvate + ions
############################################################

if [[ -d "${CHARMM36_DIR}" ]]; then
    cp -r "${CHARMM36_DIR}" .
fi

gmx_mpi pdb2gmx  -f "${INPUT_PDB}" -ignh \
    -p target_TOPOLOGY -o target \
    -ff "${FFNAME}" -water "${WATERMODEL}"

gmx_mpi editconf -f target -c -d 1.5 -bt cubic -o target_box
gmx_mpi solvate  -cp target_box -cs -p target_TOPOLOGY -o target_box_sol

gmx_mpi grompp -f "${MDP_DIR}/010_em_full_steep.mdp" \
    -c target_box_sol -p target_TOPOLOGY \
    -o targetmdrun_EM0 -r target_box_sol -maxwarn 2

echo "13" | gmx_mpi genion \
    -s targetmdrun_EM0.tpr \
    -o target_ION.gro \
    -p target_TOPOLOGY.top \
    -pname NA -nname CL -neutral -conc 0.15

if [[ -z "${INPUT_NDX}" ]]; then
    echo "q" | gmx_mpi make_ndx -f target_ION.gro -o index.ndx
    INDEX="index.ndx"
else
    INDEX="${INPUT_NDX}"
fi

############################################################
# 3) posre files
############################################################
echo "Generating posre files..."
mkdir -p posres
for posre in posre*.itp; do
    [[ -f "$posre" ]] || continue
    name=${posre%.itp}
    cp "$posre" "posres/${name}_5000.itp"
    sed -i 's/1000  1000  1000/5000  5000  5000/g' "posres/${name}_5000.itp"
    cp "$posre" "posres/${name}_1000.itp"
    cp "$posre" "posres/${name}_100.itp"
    sed -i 's/1000  1000  1000/100  100  100/g' "posres/${name}_100.itp"
    cp "$posre" "posres/${name}_10.itp"
    sed -i 's/1000  1000  1000/10  10  10/g' "posres/${name}_10.itp"
    cp "$posre" "posres/${name}_0.itp"
    sed -i 's/1000  1000  1000/0  0  0/g' "posres/${name}_0.itp"
done

update_posre() {
    local level="$1"
    for f in posres/*_"${level}".itp; do
        local base
        base=$(basename "$f")
        base=${base%_"${level}".itp}
        cp "$f" "${base}.itp"
    done
}

############################################################
# 4) Minimizacao (sem modificacoes)
############################################################

update_posre 5000
gmx_mpi grompp -f "${MDP_DIR}/005_em_pr.mdp" \
    -c target_ION.gro -r target_ION.gro \
    -p target_TOPOLOGY -o targetmdrun_EM1 -n "${INDEX}" -po mdout_EM1 -maxwarn 2
gmx_mpi mdrun -s targetmdrun_EM1 -o target_traj_EM1 -c target_BOX_EM1 \
    -g target_log_EM1 -e energy_EM1.edr -v

gmx_mpi grompp -f "${MDP_DIR}/010_em_full_steep.mdp" \
    -c target_BOX_EM1 -r target_BOX_EM1 \
    -p target_TOPOLOGY -o targetmdrun_EM2 -n "${INDEX}" -po mdout_EM2 -maxwarn 2
gmx_mpi mdrun -s targetmdrun_EM2 -o target_traj_EM2 -c target_BOX_EM2 \
    -g target_log_EM2 -e energy_EM2.edr

gmx_mpi grompp -f "${MDP_DIR}/015_em_full_cg.mdp" \
    -c target_BOX_EM2 -r target_BOX_EM2 \
    -p target_TOPOLOGY -o targetmdrun_EM3 -n "${INDEX}" -po mdout_EM3 -maxwarn 2
gmx_mpi mdrun -s targetmdrun_EM3 -o target_traj_EM3 -c target_BOX_EM3 \
    -g target_log_EM3 -e energy_EM3.edr

############################################################
# 5) Equilibracao - nova sequencia
# [A] EQ1-EQ3: NVT puro, aquecimento gradual 100->200->280K
############################################################

# EQ1: NVT 100K, posres 5000 - sem GPU (passo curto, NVT)
gmx_mpi grompp -f "${MDP_DIR}/020_md_eq_1.mdp" \
    -c target_BOX_EM3 -r target_BOX_EM3 \
    -p target_TOPOLOGY -o targetmdrun_EQ1 -n "${INDEX}" -po mdout_EQ1 -maxwarn 2
gmx_mpi mdrun -s targetmdrun_EQ1 -o target_traj_EQ1 -c target_BOX_EQ1 \
    -g target_log_EQ1 -e energy_EQ1.edr

# EQ2: NVT 200K, posres 5000
gmx_mpi grompp -f "${MDP_DIR}/025_md_eq_2.mdp" \
    -c target_BOX_EQ1 -r target_BOX_EQ1 \
    -p target_TOPOLOGY -o targetmdrun_EQ2 -n "${INDEX}" -po mdout_EQ2 -maxwarn 2
gmx_mpi mdrun -s targetmdrun_EQ2 -o target_traj_EQ2 -c target_BOX_EQ2 \
    -g target_log_EQ2 -e energy_EQ2.edr

# EQ3: NVT 280K, posres 1000
update_posre 1000  # [B] posre 1000 comeca aqui, ainda no NVT
gmx_mpi grompp -f "${MDP_DIR}/030_md_eq_3.mdp" \
    -c target_BOX_EQ2 -r target_BOX_EQ2 \
    -p target_TOPOLOGY -o targetmdrun_EQ3 -n "${INDEX}" -po mdout_EQ3 -maxwarn 2
gmx_mpi mdrun -s targetmdrun_EQ3 -o target_traj_EQ3 -c target_BOX_EQ3 \
    -g target_log_EQ3 -e energy_EQ3.edr

############################################################
# [A] EQ4-EQ6: NPT Berendsen, soltura gradual de restraints
#     [C] GPU ativado a partir daqui (passos mais longos)
############################################################

# EQ4: NPT Berendsen 280K, posres 100
update_posre 100   # [B] posre 100 no primeiro passo NPT
gmx_mpi grompp -f "${MDP_DIR}/035_md_eq_4.mdp" \
    -c target_BOX_EQ3 -r target_BOX_EQ3 \
    -p target_TOPOLOGY -o targetmdrun_EQ4 -n "${INDEX}" -po mdout_EQ4 -maxwarn 2
gmx_mpi mdrun ${PME_THREADS} -s targetmdrun_EQ4 -o target_traj_EQ4 -c target_BOX_EQ4 \
    -g target_log_EQ4 -e energy_EQ4.edr  # [C] GPU ativado

# EQ5: NPT Berendsen 280K, posres 10
update_posre 10
gmx_mpi grompp -f "${MDP_DIR}/040_md_eq_5.mdp" \
    -c target_BOX_EQ4 -r target_BOX_EQ4 \
    -p target_TOPOLOGY -o targetmdrun_EQ5 -n "${INDEX}" -po mdout_EQ5 -maxwarn 2
gmx_mpi mdrun ${PME_THREADS} -s targetmdrun_EQ5 -o target_traj_EQ5 -c target_BOX_EQ5 \
    -g target_log_EQ5 -e energy_EQ5.edr  # [C] GPU ativado

# EQ6: NPT Berendsen 280K, sem restraints
update_posre 0
gmx_mpi grompp -f "${MDP_DIR}/045_md_eq_6.mdp" \
    -c target_BOX_EQ5 -r target_BOX_EQ5 \
    -p target_TOPOLOGY -o targetmdrun_EQ6 -n "${INDEX}" -po mdout_EQ6 -maxwarn 2
gmx_mpi mdrun ${PME_THREADS} -s targetmdrun_EQ6 -o target_traj_EQ6 -c target_BOX_EQ6 \
    -g target_log_EQ6 -e energy_EQ6.edr  # [C] GPU ativado

############################################################
# [A] EQ7-EQ8: NPT Parrinello-Rahman, aquecimento 290->310K
############################################################

# EQ7: NPT PR 290K, sem restraints
gmx_mpi grompp -f "${MDP_DIR}/050_md_eq_7.mdp" \
    -c target_BOX_EQ6 -r target_BOX_EQ6 \
    -p target_TOPOLOGY -o targetmdrun_EQ7 -n "${INDEX}" -po mdout_EQ7 -maxwarn 2
gmx_mpi mdrun ${PME_THREADS} -s targetmdrun_EQ7 -o target_traj_EQ7 -c target_BOX_EQ7 \
    -g target_log_EQ7 -e energy_EQ7.edr

# EQ8: NPT PR 310K, termalizacao final
gmx_mpi grompp -f "${MDP_DIR}/055_md_eq_8.mdp" \
    -c target_BOX_EQ7 -r target_BOX_EQ7 \
    -p target_TOPOLOGY -o targetmdrun_EQ8 -n "${INDEX}" -po mdout_EQ8 -maxwarn 2
gmx_mpi mdrun ${PME_THREADS} -s targetmdrun_EQ8 -o target_traj_EQ8 -c target_BOX_EQ8 \
    -g target_log_EQ8 -e energy_EQ8.edr

############################################################
# 6) Primeiro bloco de producao: 0-5 ns
############################################################
gmx_mpi grompp -f "${MDP_DIR}/060_md_full_005ns.mdp" \
    -c target_BOX_EQ8 -r target_BOX_EQ8 \
    -p target_TOPOLOGY -o targetmdrun_MD_05ns -n "${INDEX}" -po mdout_MD_05ns.mdp -maxwarn 2
gmx_mpi mdrun ${PME_THREADS} -s targetmdrun_MD_05ns -o target_traj_MD_05ns \
    -x target_trajXTC_MD_05ns -c target_BOX_MD_05ns \
    -g target_log_MD_05ns -e energy_MD_05ns.edr

cp target_BOX_MD_05ns.gro target_BOX_MD_5ns.gro

############################################################
# 7) Loop de producao: blocos de 5 ns ate ENDTIME_NS
############################################################

CURRENT_NS=5
while (( CURRENT_NS < ENDTIME_NS )); do
    NEXT_NS=$(( CURRENT_NS + 5 ))
    echo "Running production ${CURRENT_NS} -> ${NEXT_NS} ns ..."
    cp "${MDP_DIR}/070_md_full_5ns-step.mdp" md_full.mdp

    # [D] tinit = tempo de producao acumulado
    TINIT_PS=$(( CURRENT_NS * 1000 ))
    sed -i "s/tinit               =  5000/tinit               =  ${TINIT_PS}/" md_full.mdp

    LAST_GRO="target_BOX_MD_${CURRENT_NS}ns.gro"
    if [[ ! -f "${LAST_GRO}" ]]; then
        LAST_GRO="target_BOX_MD_${CURRENT_NS}.gro"
    fi

    gmx_mpi grompp -f md_full.mdp \
        -c "${LAST_GRO}" -r "${LAST_GRO}" \
        -p target_TOPOLOGY -o "targetmdrun_MD_${NEXT_NS}ns" \
        -po "mdout_MD_${NEXT_NS}ns.mdp" -n "${INDEX}" -maxwarn 2

    gmx_mpi mdrun ${PME_THREADS} \
        -s "targetmdrun_MD_${NEXT_NS}ns" \
        -o "target_traj_MD_${NEXT_NS}ns" \
        -x "target_trajXTC_MD_${NEXT_NS}ns" \
        -c "target_BOX_MD_${NEXT_NS}ns" \
        -g "target_log_MD_${NEXT_NS}ns" \
        -e "energy_MD_${NEXT_NS}ns.edr"

    CURRENT_NS=${NEXT_NS}
done

rm -f md_full.mdp

ENDTIME=$(date "+%H:%M:%S-%D")
echo "Run ended at ${ENDTIME}"

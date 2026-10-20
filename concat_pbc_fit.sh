#!/bin/bash
# Uso: ./concat_pbc_fit.sh [NOME]
# Rodar na pasta DM_0. Usa target_trajXTC_MD_*ns.xtc (05ns ... 500ns).
#
# Saidas:
#   trjCAT_NOME_05ns-500ns.xtc                           (1) trajetorias concatenadas
#   trjCONV_NOME_05-500ns_pbc-fit-WithWater.xtc          (2)+(3) pbc nojump + fit rot+trans, com agua
#   trjCONV_NOME_05-500ns_pbc-fit-WithOUTWater.xtc       (4) mesmo, sem agua (ligante e ions ficam)
#
# Opcionais: GMX=gmx_mpi  REF=targetmdrun_EQ1.tpr  FIT_GROUP=Protein  NOWATER_GROUP=non-Water
set -euo pipefail

GMX="${GMX:-gmx_mpi}"
REF="${REF:-targetmdrun_EQ1.tpr}"
FIT_GROUP="${FIT_GROUP:-Protein}"
NOWATER_GROUP="${NOWATER_GROUP:-non-Water}"
name="${1:-bes_ligand}"

mapfile -t catMDs < <(ls target_trajXTC_MD_*ns.xtc | sort -V)
[ ${#catMDs[@]} -gt 0 ] || { echo "Nenhum target_trajXTC_MD_*ns.xtc aqui."; exit 1; }
[ -s "$REF" ] || { echo "Nao achei $REF"; exit 1; }
first=$(sed -E 's/.*_MD_([0-9]+)ns\.xtc/\1/' <<<"${catMDs[0]}")
last=$(sed -E 's/.*_MD_([0-9]+)ns\.xtc/\1/' <<<"${catMDs[-1]}")
echo "${#catMDs[@]} segmentos: ${first}ns -> ${last}ns | referencia: $REF"

catFULL=trjCAT_"$name"_"$first"ns-"$last"ns.xtc
pbc=trjCONV_"$name"_"$first"-"$last"ns_pbc.xtc
fitW=trjCONV_"$name"_"$first"-"$last"ns_pbc-fit-WithWater.xtc
fitNW=trjCONV_"$name"_"$first"-"$last"ns_pbc-fit-WithOUTWater.xtc

# 1) concatenar
$GMX trjcat -f "${catMDs[@]}" -o "$catFULL"
echo "Finished file $catFULL"

# 2) remover saltos de PBC (saida: System)
printf 'System\n' | $GMX trjconv -f "$catFULL" -s "$REF" -o "$pbc" -pbc nojump

# 3) ajuste rot+trans (fit na proteina, saida: System = com agua)
printf '%s\nSystem\n' "$FIT_GROUP" | $GMX trjconv -f "$pbc" -s "$REF" -o "$fitW" -fit rot+trans
rm "$pbc"

# 4) versao sem agua (proteina + ligante + ions)
printf '%s\n' "$NOWATER_GROUP" | $GMX trjconv -f "$fitW" -s "$REF" -o "$fitNW"

echo "Pronto:"
ls -lh "$catFULL" "$fitW" "$fitNW"

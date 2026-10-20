#!/usr/bin/env bash
set -euo pipefail

export PATH="/usr/local/gmx_plumed/bin:$PATH"

############################################################
# GROMACS + CHARMM36/CGenFF MD pipeline
# Protein + DockThor 2,4-D ligand
#
# System:
#   Protein: CHARMM36
#   Ligand:  CGenFF / DockThor 2,4-D
#   Water:   CHARMM36 TIP3P topology
#
# Important:
#   spc216.gro is used as the 3-point water coordinate template.
#   GROMACS supports spc216.gro for TIP3P-like 3-point water models.
#
# Pipeline:
#
#   1. Validate inputs
#   2. Clean generated files
#   3. Convert ligand PDB -> GRO
#   4. Generate protein topology with pdb2gmx
#   5. Extract complete CGenFF ligand moleculetype
#   6. Integrate ligand topology into protein topology
#   7. Combine protein + ligand coordinates
#   8. Define simulation box
#   9. Solvate
#  10. Create pre-ionization TPR
#  11. Add 0.15 M NaCl + neutralize
#  12. Create index
#  13. Generate protein position restraints
#  14. Energy minimization
#  15. Equilibration
#  16. Production to ENDTIME_NS
#
# The only grompp call using -maxwarn is the PRE-IONIZATION
# grompp. At that point the system is intentionally charged:
#
#   protein = +4
#   ligand  = -1
#   ----------------
#   system  = +3
#
# This warning is expected because genion has not yet added
# counterions.
############################################################

set -o errtrace

############################################################
# USER SETTINGS
############################################################

ENDTIME_NS=500

MDP_DIR="/home/lorenzo/DM_iGEM/MDPs_novos"

CHARMM36_DIR="charmm36.ff"
WATERMODEL="tip3p"
FFNAME="charmm36"

# GPU settings used for the later MD stages
PME_THREADS="-nb gpu -bonded gpu -pme gpu -pin on -nstlist 400"

PROTEIN_PDB="protein_a8119c2bf1_prep_fixed.pdb"
LIGAND_PDB="bestranking_gmx.pdb"
LIGAND_TOP="bestranking_gmx.top"

PROTEIN_GRO="target_protein.gro"
COMPLEX_GRO="target_complex.gro"
LIGAND_GRO="ligand.gro"

TOPOLOGY="target_TOPOLOGY.top"
LIGAND_ITP="bes_ligand.itp"
INDEX="index.ndx"

WATER_GRO="${CHARMM36_DIR}/spc216.gro"

############################################################
# FUNCTIONS
############################################################

log() {
    echo
    echo "============================================================"
    echo "$*"
    echo "============================================================"
}

require_file() {
    local f="$1"

    if [[ ! -f "$f" ]]; then
        echo "ERROR: required file not found: $f" >&2
        exit 1
    fi
}

require_dir() {
    local d="$1"

    if [[ ! -d "$d" ]]; then
        echo "ERROR: required directory not found: $d" >&2
        exit 1
    fi
}

run_grompp() {
    gmx_mpi grompp "$@"
}

############################################################
# START TIME
############################################################

STARTTIME=$(date "+%H:%M:%S-%D")

echo "Run started at ${STARTTIME}"

############################################################
# 0) VALIDATE INPUTS
############################################################

log "0) Validating input files"

require_file "${PROTEIN_PDB}"
require_file "${LIGAND_PDB}"
require_file "${LIGAND_TOP}"

require_dir "${CHARMM36_DIR}"

require_file "${CHARMM36_DIR}/forcefield.itp"
require_file "${CHARMM36_DIR}/ffnonbonded.itp"
require_file "${CHARMM36_DIR}/ffbonded.itp"
require_file "${CHARMM36_DIR}/tip3p.itp"
require_file "${CHARMM36_DIR}/ions.itp"
require_file "${CHARMM36_DIR}/bes_ffbonded.itp"
require_file "${WATER_GRO}"

for mdp in \
    005_em_pr.mdp \
    010_em_full_steep.mdp \
    015_em_full_cg.mdp \
    020_md_eq_1.mdp \
    025_md_eq_2.mdp \
    030_md_eq_3.mdp \
    035_md_eq_4.mdp \
    040_md_eq_5.mdp \
    045_md_eq_6.mdp \
    050_md_eq_7.mdp \
    055_md_eq_8.mdp \
    060_md_full_005ns.mdp \
    070_md_full_5ns-step.mdp
do
    require_file "${MDP_DIR}/${mdp}"
done

echo
echo "Protein:       ${PROTEIN_PDB}"
echo "Ligand PDB:    ${LIGAND_PDB}"
echo "Ligand TOP:    ${LIGAND_TOP}"
echo "CHARMM36 dir:  ${CHARMM36_DIR}"
echo "Water coords:  ${WATER_GRO}"
echo "MDP directory: ${MDP_DIR}"
echo "End time:      ${ENDTIME_NS} ns"

############################################################
# 1) CLEAN GENERATED FILES
############################################################

log "1) Cleaning generated files"

rm -f \
    target \
    target.gro \
    target_protein.gro \
    target_complex.gro \
    ligand.gro \
    target_box.gro \
    target_box_sol.gro \
    target_ION.gro \
    target_TOPOLOGY.top \
    bes_ligand.itp \
    index.ndx \
    posre.itp \
    md_full.mdp \
    mdout.mdp \
    mdout_*.mdp \
    targetmdrun_*.tpr \
    targetmdrun_* \
    target_traj_* \
    target_trajXTC_* \
    target_log_* \
    target_BOX_* \
    energy_*.edr \
    mdout_* \
    \#posre.itp.* \
    \#target_TOPOLOGY.top.*

rm -rf posres

############################################################
# 2) CONVERT LIGAND PDB -> GRO
############################################################

log "2) Converting CGenFF ligand PDB to GRO"

gmx_mpi editconf \
    -f "${LIGAND_PDB}" \
    -o "${LIGAND_GRO}"

############################################################
# 3) BUILD PROTEIN TOPOLOGY
############################################################

log "3) Generating protein topology with pdb2gmx"

gmx_mpi pdb2gmx \
    -f "${PROTEIN_PDB}" \
    -ignh \
    -p "${TOPOLOGY}" \
    -o "${PROTEIN_GRO}" \
    -ff "${FFNAME}" \
    -water "${WATERMODEL}"

############################################################
# 4) EXTRACT COMPLETE CGenFF LIGAND MOLECULETYPE
############################################################
#
# Preserve the complete CGenFF [ moleculetype ] section:
#
#   [ moleculetype ]
#   [ atoms ]
#   [ bonds ]
#   [ angles ]
#   [ dihedrals ]
#   [ impropers ]
#   [ pairs ]
#   [ exclusions ]
#   [ virtual_sites3 ]
#   etc.
#
# Stop before the system-level [ system ] section.
#
# IMPORTANT:
# The original CGenFF topology may contain includes for:
#   - posre.itp
#   - tip3p.itp
#   - ions.itp
#
# These must NOT be included inside the ligand topology.
# They are system-level/protein-level includes and are handled
# separately by target_TOPOLOGY.top.
############################################################

log "4) Extracting complete CGenFF ligand topology"

python3 - "${LIGAND_TOP}" "${LIGAND_ITP}" <<'PYCODE'
import sys
import re
from pathlib import Path

src = Path(sys.argv[1])
out = Path(sys.argv[2])

lines = src.read_text().splitlines()

# ------------------------------------------------------------
# Find ligand [ moleculetype ]
# ------------------------------------------------------------

start = None

for i, line in enumerate(lines):
    if line.strip().lower() == "[ moleculetype ]":
        start = i
        break

if start is None:
    raise SystemExit(
        "ERROR: [ moleculetype ] was not found in the CGenFF topology."
    )

# ------------------------------------------------------------
# Find the first system-level [ system ] after the ligand
# ------------------------------------------------------------

system_idx = None

for i in range(start + 1, len(lines)):
    if lines[i].strip().lower() == "[ system ]":
        system_idx = i
        break

if system_idx is None:
    raise SystemExit(
        "ERROR: [ system ] was not found after ligand [ moleculetype ]."
    )

# Keep everything belonging to the ligand moleculetype.
section = lines[start:system_idx]

# ------------------------------------------------------------
# Remove system-level includes
# ------------------------------------------------------------
#
# The CGenFF-generated topology may contain:
#
#   #include "posre.itp"
#   #include "./charmm36.ff/tip3p.itp"
#   #include "./charmm36.ff/ions.itp"
#
# None of these belong inside bes_ligand.itp.
#
# The main target topology handles:
#   - protein position restraints
#   - TIP3P water
#   - ions
#
# The ligand itself remains completely unrestrained.
# ------------------------------------------------------------

filtered = []

removed_includes = []

for line in section:

    stripped = line.strip()

    if stripped.lower().startswith("#include"):

        lower = stripped.lower()

        if (
            "posre" in lower
            or "tip3p.itp" in lower
            or "ions.itp" in lower
        ):
            removed_includes.append(stripped)
            continue

    filtered.append(line)

section = filtered

# ------------------------------------------------------------
# Convert to text
# ------------------------------------------------------------

text = "\n".join(section)

# ------------------------------------------------------------
# Basic structural checks
# ------------------------------------------------------------

required_sections = [
    "[ moleculetype ]",
    "[ atoms ]",
    "[ bonds ]",
]

text_lower = text.lower()

for required in required_sections:
    if required.lower() not in text_lower:
        raise SystemExit(
            f"ERROR: extracted ligand topology does not contain {required}."
        )

# ------------------------------------------------------------
# Make sure no forbidden system-level includes survived
# ------------------------------------------------------------

for line in section:
    lower = line.strip().lower()

    if lower.startswith("#include"):
        if (
            "posre" in lower
            or "tip3p.itp" in lower
            or "ions.itp" in lower
        ):
            raise SystemExit(
                f"ERROR: forbidden system-level include survived extraction:\n"
                f"    {line}"
            )

# ------------------------------------------------------------
# Validate ligand moleculetype name
# ------------------------------------------------------------

m = re.search(
    r"\[\s*moleculetype\s*\](.*?)(?=\n\s*\[|\Z)",
    text,
    re.I | re.S,
)

if not m:
    raise SystemExit(
        "ERROR: could not parse ligand [ moleculetype ] section."
    )

molecule_lines = [
    line.strip()
    for line in m.group(1).splitlines()
    if line.strip() and not line.strip().startswith(";")
]

if not molecule_lines:
    raise SystemExit(
        "ERROR: ligand [ moleculetype ] has no molecule name."
    )

ligand_moltype = molecule_lines[0].split()[0]

# ------------------------------------------------------------
# Write output
# ------------------------------------------------------------

header = (
    "; ============================================================\n"
    "; CGenFF ligand topology extracted from bestranking_gmx.top\n"
    "; ============================================================\n"
    "; The CHARMM36 force field is included by target_TOPOLOGY.top.\n"
    "; CGenFF-specific bonded parameters are included separately.\n"
    "; System-level includes were intentionally removed:\n"
    ";   - posre.itp\n"
    ";   - tip3p.itp\n"
    ";   - ions.itp\n"
    "; The ligand is intentionally not position restrained.\n"
    ";\n"
)

out.write_text(
    header
    + text
    + "\n"
)

print(f"Ligand moleculetype: {ligand_moltype}")
print(f"Wrote: {out}")

if removed_includes:
    print("Removed system-level includes:")
    for include in removed_includes:
        print(f"  {include}")
else:
    print("No system-level includes needed to be removed.")
PYCODE

############################################################
# 5) VALIDATE LIGAND TOPOLOGY
############################################################

log "5) Validating ligand topology"

python3 - "${LIGAND_ITP}" <<'PYCODE'
import sys
import re
from pathlib import Path

path = Path(sys.argv[1])
text = path.read_text()

# Locate atoms section.
m = re.search(
    r"\[\s*atoms\s*\](.*?)(?=\n\s*\[|\Z)",
    text,
    re.I | re.S,
)

if not m:
    raise SystemExit(
        "ERROR: could not find ligand [ atoms ] section."
    )

rows = []

for line in m.group(1).splitlines():

    line = line.split(";", 1)[0].strip()

    if not line:
        continue

    if line.startswith("#"):
        continue

    fields = line.split()

    if fields and fields[0].isdigit():
        rows.append(fields)

expected_atoms = 20

if len(rows) != expected_atoms:
    raise SystemExit(
        f"ERROR: expected {expected_atoms} ligand atoms "
        f"including virtual sites; found {len(rows)}."
    )

# GROMACS [ atoms ]:
#
# nr type resnr residue atom cgnr charge mass
#
charge = sum(float(row[6]) for row in rows)

print(f"Ligand topology atoms:  {len(rows)}")
print(f"Ligand topology charge: {charge:.6f}")

if abs(charge + 1.0) > 1e-5:
    raise SystemExit(
        "ERROR: ligand topology total charge is not -1.0."
    )

print("Ligand topology validation: OK")
PYCODE

############################################################
# 6) INTEGRATE LIGAND INTO SYSTEM TOPOLOGY
############################################################

log "6) Integrating ligand topology into system topology"

python3 - "${TOPOLOGY}" "${LIGAND_ITP}" "${CHARMM36_DIR}" <<'PYCODE'
import sys
from pathlib import Path

top_path = Path(sys.argv[1])
ligand_itp = Path(sys.argv[2])
ffdir = sys.argv[3]

lines = top_path.read_text().splitlines()

############################################################
# Add CGenFF-specific bonded parameters
############################################################

bonded_include = f'#include "{ffdir}/bes_ffbonded.itp"'

if not any("bes_ffbonded.itp" in line for line in lines):

    insert_at = None

    for i, line in enumerate(lines):

        if (
            "forcefield.itp" in line
            and line.lstrip().startswith("#include")
        ):
            insert_at = i + 1
            break

    if insert_at is None:
        raise SystemExit(
            "ERROR: could not find forcefield.itp include."
        )

    lines.insert(insert_at, bonded_include)

############################################################
# Add ligand moleculetype include
############################################################

if any(ligand_itp.name in line for line in lines):
    raise SystemExit(
        "ERROR: ligand topology is already included."
    )

system_idx = None

for i, line in enumerate(lines):

    if line.strip().lower() == "[ system ]":
        system_idx = i
        break

if system_idx is None:
    raise SystemExit(
        "ERROR: [ system ] not found."
    )

lines.insert(
    system_idx,
    f'#include "{ligand_itp.name}"'
)

############################################################
# Replace system name
############################################################

for i, line in enumerate(lines):

    if line.strip().lower() == "[ system ]":

        if i + 1 < len(lines):
            lines[i + 1] = "Protein + 2,4-D (DockThor/CGenFF)"

        break

############################################################
# Add ligand to [ molecules ]
############################################################

molecules_idx = None

for i, line in enumerate(lines):

    if line.strip().lower() == "[ molecules ]":
        molecules_idx = i
        break

if molecules_idx is None:
    raise SystemExit(
        "ERROR: [ molecules ] not found."
    )

# Make sure we don't already have Other.
for line in lines[molecules_idx + 1:]:

    stripped = line.strip()

    if stripped.startswith("["):
        break

    if stripped.split()[:1] == ["Other"]:
        raise SystemExit(
            "ERROR: Other is already present in [ molecules ]."
        )

# Insert after existing molecule entries and before next section.
insert_idx = len(lines)

for i in range(molecules_idx + 1, len(lines)):

    if lines[i].strip().startswith("["):
        insert_idx = i
        break

lines.insert(
    insert_idx,
    "Other                 1"
)

top_path.write_text(
    "\n".join(lines) + "\n"
)
PYCODE

echo
echo "Topology includes:"
grep -n \
    "bes_ffbonded.itp\|bes_ligand.itp\|Other" \
    "${TOPOLOGY}"

############################################################
# 7) COMBINE PROTEIN + LIGAND COORDINATES
############################################################

log "7) Combining protein and DockThor ligand coordinates"

python3 - \
    "${PROTEIN_GRO}" \
    "${LIGAND_GRO}" \
    "${COMPLEX_GRO}" <<'PYCODE'

import sys
from pathlib import Path

protein_path = Path(sys.argv[1])
ligand_path = Path(sys.argv[2])
output_path = Path(sys.argv[3])


def read_gro(path):

    lines = path.read_text().splitlines()

    if len(lines) < 3:
        raise SystemExit(
            f"ERROR: invalid GRO file: {path}"
        )

    try:
        natoms = int(lines[1].strip())
    except ValueError:
        raise SystemExit(
            f"ERROR: invalid atom count in {path}"
        )

    atom_lines = lines[2:2 + natoms]

    if len(atom_lines) != natoms:
        raise SystemExit(
            f"ERROR: atom count mismatch in {path}: "
            f"expected {natoms}, found {len(atom_lines)}"
        )

    box_line = lines[2 + natoms]

    return atom_lines, box_line


protein_atoms, protein_box = read_gro(protein_path)
ligand_atoms, _ = read_gro(ligand_path)

############################################################
# Expected atom counts
############################################################

if len(protein_atoms) != 5288:
    raise SystemExit(
        f"ERROR: expected 5288 protein atoms after pdb2gmx; "
        f"found {len(protein_atoms)}."
    )

if len(ligand_atoms) != 20:
    raise SystemExit(
        f"ERROR: expected 20 ligand atoms including virtual sites; "
        f"found {len(ligand_atoms)}."
    )

############################################################
# Renumber ligand atoms
############################################################

new_ligand_atoms = []

first_ligand_atom = len(protein_atoms) + 1

for i, line in enumerate(ligand_atoms, start=0):

    if len(line) < 44:
        raise SystemExit(
            "ERROR: malformed ligand GRO atom line."
        )

    resname = line[5:10]
    atomname = line[10:15]
    coords = line[20:44]

    velocities = line[44:] if len(line) > 44 else ""

    atom_number = first_ligand_atom + i

    new_line = (
        f"{1:5d}"
        f"{resname:>5}"
        f"{atomname:>5}"
        f"{atom_number:5d}"
        f"{coords}"
        f"{velocities}"
    )

    new_ligand_atoms.append(new_line)

############################################################
# Build combined structure
############################################################

combined = protein_atoms + new_ligand_atoms

if len(combined) != 5308:
    raise SystemExit(
        f"ERROR: expected 5308 total atoms; "
        f"found {len(combined)}."
    )

output_path.write_text(
    "Protein + DockThor 2,4-D\n"
    f"{len(combined):5d}\n"
    + "\n".join(combined)
    + "\n"
    + protein_box
    + "\n"
)

print(f"Protein atoms: {len(protein_atoms)}")
print(f"Ligand atoms:  {len(ligand_atoms)}")
print(f"Complex atoms: {len(combined)}")
PYCODE

############################################################
# 8) DEFINE BOX
############################################################

log "8) Defining simulation box"

gmx_mpi editconf \
    -f "${COMPLEX_GRO}" \
    -c \
    -d 1.5 \
    -bt cubic \
    -o target_box.gro

############################################################
# 9) SOLVATE
############################################################

log "9) Solvating system"

gmx_mpi solvate \
    -cp target_box.gro \
    -cs "${WATER_GRO}" \
    -p "${TOPOLOGY}" \
    -o target_box_sol.gro

############################################################
# 10) PREPARE TPR FOR GENION
############################################################
#
# At this point:
#
#   Protein = +4
#   Ligand  = -1
#   ----------------
#   Total   = +3
#
# This is intentional.
#
# grompp otherwise stops because PME/Ewald requires a neutral
# system. We allow exactly this one known warning.
#
# No other grompp call below uses -maxwarn.
############################################################

log "10) Preparing solvated system for ion addition"

gmx_mpi grompp \
    -f "${MDP_DIR}/010_em_full_steep.mdp" \
    -c target_box_sol.gro \
    -p "${TOPOLOGY}" \
    -o targetmdrun_EM0.tpr \
    -r target_box_sol.gro \
    -maxwarn 1

############################################################
# 11) ADD NaCl + NEUTRALIZE
############################################################

log "11) Adding 0.15 M NaCl and neutralizing system"

echo "SOL" | gmx_mpi genion \
    -s targetmdrun_EM0.tpr \
    -o target_ION.gro \
    -p "${TOPOLOGY}" \
    -pname NA \
    -nname CL \
    -neutral \
    -conc 0.15

############################################################
# VERIFY IONIZATION OUTPUT
############################################################

log "Checking ionized topology"

grep -A20 \
    '^\[ molecules \]' \
    "${TOPOLOGY}" \
    || true

############################################################
# 12) CREATE FRESH INDEX
############################################################

log "12) Creating fresh index"

echo -e "r bes\nq" | \
    gmx_mpi make_ndx \
        -f target_ION.gro \
        -o "${INDEX}"

############################################################
# 13) GENERATE POSITION RESTRAINT VARIANTS
############################################################

log "13) Generating protein position-restraint variants"

if [[ ! -f posre.itp ]]; then
    echo "ERROR: pdb2gmx did not create posre.itp." >&2
    exit 1
fi

mkdir -p posres

cp posre.itp posres/posre_5000.itp

sed -i \
    's/1000[[:space:]]\+1000[[:space:]]\+1000/5000  5000  5000/g' \
    posres/posre_5000.itp

cp posre.itp posres/posre_1000.itp

cp posre.itp posres/posre_100.itp

sed -i \
    's/1000[[:space:]]\+1000[[:space:]]\+1000/100  100  100/g' \
    posres/posre_100.itp

cp posre.itp posres/posre_10.itp

sed -i \
    's/1000[[:space:]]\+1000[[:space:]]\+1000/10  10  10/g' \
    posres/posre_10.itp

cp posre.itp posres/posre_0.itp

sed -i \
    's/1000[[:space:]]\+1000[[:space:]]\+1000/0  0  0/g' \
    posres/posre_0.itp


update_posre() {

    local level="$1"

    cp \
        "posres/posre_${level}.itp" \
        posre.itp
}

############################################################
# 14) ENERGY MINIMIZATION
############################################################

log "14.1) EM1: restrained minimization"

update_posre 5000

run_grompp \
    -f "${MDP_DIR}/005_em_pr.mdp" \
    -c target_ION.gro \
    -r target_ION.gro \
    -p "${TOPOLOGY}" \
    -o targetmdrun_EM1.tpr \
    -n "${INDEX}" \
    -po mdout_EM1.mdp \
    -maxwarn 2

gmx_mpi mdrun \
    -s targetmdrun_EM1.tpr \
    -o target_traj_EM1 \
    -c target_BOX_EM1 \
    -g target_log_EM1 \
    -e energy_EM1.edr \
    -v

############################################################

log "14.2) EM2: steepest descent"

run_grompp \
    -f "${MDP_DIR}/010_em_full_steep.mdp" \
    -c target_BOX_EM1.gro \
    -r target_BOX_EM1.gro \
    -p "${TOPOLOGY}" \
    -o targetmdrun_EM2.tpr \
    -n "${INDEX}" \
    -po mdout_EM2.mdp \
    -maxwarn 2

gmx_mpi mdrun \
    -s targetmdrun_EM2.tpr \
    -o target_traj_EM2 \
    -c target_BOX_EM2 \
    -g target_log_EM2 \
    -e energy_EM2.edr

############################################################

log "14.3) EM3: conjugate gradient"

run_grompp \
    -f "${MDP_DIR}/015_em_full_cg.mdp" \
    -c target_BOX_EM2.gro \
    -r target_BOX_EM2.gro \
    -p "${TOPOLOGY}" \
    -o targetmdrun_EM3.tpr \
    -n "${INDEX}" \
    -po mdout_EM3.mdp \
    -maxwarn 2

gmx_mpi mdrun \
    -s targetmdrun_EM3.tpr \
    -o target_traj_EM3 \
    -c target_BOX_EM3 \
    -g target_log_EM3 \
    -e energy_EM3.edr

############################################################
# 15) EQUILIBRATION
############################################################

log "15.1) EQ1: NVT 100 K, protein posres 5000"

update_posre 5000

run_grompp \
    -f "${MDP_DIR}/020_md_eq_1.mdp" \
    -c target_BOX_EM3.gro \
    -r target_BOX_EM3.gro \
    -p "${TOPOLOGY}" \
    -o targetmdrun_EQ1.tpr \
    -n "${INDEX}" \
    -po mdout_EQ1.mdp \
    -maxwarn 2

gmx_mpi mdrun \
    -s targetmdrun_EQ1.tpr \
    -o target_traj_EQ1 \
    -c target_BOX_EQ1 \
    -g target_log_EQ1 \
    -e energy_EQ1.edr

############################################################

log "15.2) EQ2: NVT 200 K, protein posres 5000"

run_grompp \
    -f "${MDP_DIR}/025_md_eq_2.mdp" \
    -c target_BOX_EQ1.gro \
    -r target_BOX_EQ1.gro \
    -p "${TOPOLOGY}" \
    -o targetmdrun_EQ2.tpr \
    -n "${INDEX}" \
    -po mdout_EQ2.mdp \
    -maxwarn 2

gmx_mpi mdrun \
    -s targetmdrun_EQ2.tpr \
    -o target_traj_EQ2 \
    -c target_BOX_EQ2 \
    -g target_log_EQ2 \
    -e energy_EQ2.edr

############################################################

log "15.3) EQ3: NVT 280 K, protein posres 1000"

update_posre 1000

run_grompp \
    -f "${MDP_DIR}/030_md_eq_3.mdp" \
    -c target_BOX_EQ2.gro \
    -r target_BOX_EQ2.gro \
    -p "${TOPOLOGY}" \
    -o targetmdrun_EQ3.tpr \
    -n "${INDEX}" \
    -po mdoutEQ3.mdp \
    -maxwarn 2

gmx_mpi mdrun \
    -s targetmdrun_EQ3.tpr \
    -o target_traj_EQ3 \
    -c target_BOX_EQ3 \
    -g target_log_EQ3 \
    -e energy_EQ3.edr

############################################################

log "15.4) EQ4: NPT Berendsen 280 K, protein posres 100"

update_posre 100

run_grompp \
    -f "${MDP_DIR}/035_md_eq_4.mdp" \
    -c target_BOX_EQ3.gro \
    -r target_BOX_EQ3.gro \
    -p "${TOPOLOGY}" \
    -o targetmdrun_EQ4.tpr \
    -n "${INDEX}" \
    -po mdout_EQ4.mdp \
    -maxwarn 2

gmx_mpi mdrun \
    ${PME_THREADS} \
    -s targetmdrun_EQ4.tpr \
    -o target_traj_EQ4 \
    -c target_BOX_EQ4 \
    -g target_log_EQ4 \
    -e energy_EQ4.edr

############################################################

log "15.5) EQ5: NPT Berendsen 280 K, protein posres 10"

update_posre 10

run_grompp \
    -f "${MDP_DIR}/040_md_eq_5.mdp" \
    -c target_BOX_EQ4.gro \
    -r target_BOX_EQ4.gro \
    -p "${TOPOLOGY}" \
    -o targetmdrun_EQ5.tpr \
    -n "${INDEX}" \
    -po mdout_EQ5.mdp \
    -maxwarn 2

gmx_mpi mdrun \
    ${PME_THREADS} \
    -s targetmdrun_EQ5.tpr \
    -o target_traj_EQ5 \
    -c target_BOX_EQ5 \
    -g target_log_EQ5 \
    -e energy_EQ5.edr

############################################################

log "15.6) EQ6: NPT Berendsen 280 K, no protein restraints"

update_posre 0

run_grompp \
    -f "${MDP_DIR}/045_md_eq_6.mdp" \
    -c target_BOX_EQ5.gro \
    -r target_BOX_EQ5.gro \
    -p "${TOPOLOGY}" \
    -o targetmdrun_EQ6.tpr \
    -n "${INDEX}" \
    -po mdout_EQ6.mdp \
    -maxwarn 2

gmx_mpi mdrun \
    ${PME_THREADS} \
    -s targetmdrun_EQ6.tpr \
    -o target_traj_EQ6 \
    -c target_BOX_EQ6 \
    -g target_log_EQ6 \
    -e energy_EQ6.edr

############################################################

log "15.7) EQ7: NPT Parrinello-Rahman 290 K"

run_grompp \
    -f "${MDP_DIR}/050_md_eq_7.mdp" \
    -c target_BOX_EQ6.gro \
    -r target_BOX_EQ6.gro \
    -p "${TOPOLOGY}" \
    -o targetmdrun_EQ7.tpr \
    -n "${INDEX}" \
    -po mdout_EQ7.mdp \
    -maxwarn 2

gmx_mpi mdrun \
    ${PME_THREADS} \
    -s targetmdrun_EQ7.tpr \
    -o target_traj_EQ7 \
    -c target_BOX_EQ7 \
    -g target_log_EQ7 \
    -e energy_EQ7.edr

############################################################

log "15.8) EQ8: NPT Parrinello-Rahman 310 K"

run_grompp \
    -f "${MDP_DIR}/055_md_eq_8.mdp" \
    -c target_BOX_EQ7.gro \
    -r target_BOX_EQ7.gro \
    -p "${TOPOLOGY}" \
    -o targetmdrun_EQ8.tpr \
    -n "${INDEX}" \
    -po mdout_EQ8.mdp \
    -maxwarn 2

gmx_mpi mdrun \
    ${PME_THREADS} \
    -s targetmdrun_EQ8.tpr \
    -o target_traj_EQ8 \
    -c target_BOX_EQ8 \
    -g target_log_EQ8 \
    -e energy_EQ8.edr

############################################################
# 16) PRODUCTION: 0 -> 5 ns
############################################################

log "16) Production: 0 -> 5 ns"

run_grompp \
    -f "${MDP_DIR}/060_md_full_005ns.mdp" \
    -c target_BOX_EQ8.gro \
    -r target_BOX_EQ8.gro \
    -p "${TOPOLOGY}" \
    -o targetmdrun_MD_05ns.tpr \
    -n "${INDEX}" \
    -po mdout_MD_05ns.mdp \
    -maxwarn 2

gmx_mpi mdrun \
    ${PME_THREADS} \
    -s targetmdrun_MD_05ns.tpr \
    -o target_traj_MD_05ns \
    -x target_trajXTC_MD_05ns \
    -c target_BOX_MD_05ns \
    -g target_log_MD_05ns \
    -e energy_MD_05ns.edr

cp \
    target_BOX_MD_05ns.gro \
    target_BOX_MD_5ns.gro

############################################################
# 17) PRODUCTION: 5 -> ENDTIME_NS
############################################################

log "17) Production loop"

CURRENT_NS=5

while (( CURRENT_NS < ENDTIME_NS )); do

    NEXT_NS=$((CURRENT_NS + 5))

    if (( NEXT_NS > ENDTIME_NS )); then
        NEXT_NS=${ENDTIME_NS}
    fi

    echo
    echo "Running production ${CURRENT_NS} -> ${NEXT_NS} ns ..."

    ########################################################
    # Prepare production MDP
    ########################################################

    cp \
        "${MDP_DIR}/070_md_full_5ns-step.mdp" \
        md_full.mdp

    TINIT_PS=$((CURRENT_NS * 1000))

    sed -i -E \
        "s/^([[:space:]]*tinit[[:space:]]*=[[:space:]]*)[0-9.+-]+/\1${TINIT_PS}/" \
        md_full.mdp

    ########################################################
    # Find previous structure
    ########################################################

    LAST_GRO="target_BOX_MD_${CURRENT_NS}ns.gro"

    if [[ ! -f "${LAST_GRO}" ]]; then
        LAST_GRO="target_BOX_MD_${CURRENT_NS}.gro"
    fi

    if [[ ! -f "${LAST_GRO}" ]]; then
        echo
        echo "ERROR: production input structure not found:"
        echo "       ${LAST_GRO}"
        exit 1
    fi

    ########################################################
    # Prepare TPR
    ########################################################

    run_grompp \
        -f md_full.mdp \
        -c "${LAST_GRO}" \
        -r "${LAST_GRO}" \
        -p "${TOPOLOGY}" \
        -o "targetmdrun_MD_${NEXT_NS}ns.tpr" \
        -po "mdout_MD_${NEXT_NS}ns.mdp" \
        -n "${INDEX}" \
        -maxwarn 2

    ########################################################
    # Run production
    ########################################################

    gmx_mpi mdrun \
        ${PME_THREADS} \
        -s "targetmdrun_MD_${NEXT_NS}ns.tpr" \
        -o "target_traj_MD_${NEXT_NS}ns" \
        -x "target_trajXTC_MD_${NEXT_NS}ns" \
        -c "target_BOX_MD_${NEXT_NS}ns" \
        -g "target_log_MD_${NEXT_NS}ns" \
        -e "energy_MD_${NEXT_NS}ns.edr"

    CURRENT_NS=${NEXT_NS}
done

############################################################
# 18) CLEAN TEMPORARY PRODUCTION MDP
############################################################

rm -f md_full.mdp

############################################################
# FINISHED
############################################################

ENDTIME=$(date "+%H:%M:%S-%D")

echo
echo "============================================================"
echo "MD pipeline completed successfully."
echo
echo "Started: ${STARTTIME}"
echo "Ended:   ${ENDTIME}"
echo "Total production target: ${ENDTIME_NS} ns"
echo "============================================================"


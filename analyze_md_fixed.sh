#!/usr/bin/env bash
# Protein–ligand MD analysis pipeline for GROMACS.
# Run from the simulation directory containing targetmdrun_MD_500ns.tpr,
# target_trajXTC_MD_500ns.xtc, and energy_MD_500ns.edr.
#
# Usage:
#   chmod +x analyze_md.sh
#   ./analyze_md.sh
#
# Optional environment overrides:
#   GMX=gmx_mpi LIGAND_RESNAME=Other REF_TPR=targetmdrun_MD_500ns.tpr \
#   TRAJ=target_trajXTC_MD_500ns.xtc ENERGY=energy_MD_500ns.edr ./analyze_md.sh

set -Eeuo pipefail
IFS=$'\n\t'

GMX="${GMX:-gmx_mpi}"
OUT="${OUT:-analysis}"
REF_TPR="${REF_TPR:-targetmdrun_MD_500ns.tpr}"
TRAJ="${TRAJ:-target_trajXTC_MD_500ns.xtc}"
ENERGY="${ENERGY:-energy_MD_500ns.edr}"
STRUCTURE="${STRUCTURE:-target_BOX_MD_500ns.gro}"
LIGAND_RESNAME="${LIGAND_RESNAME:-Other}"
CUTOFF_NM="${CUTOFF_NM:-0.4}"

mkdir -p "$OUT"/{trajectories,structures,index,rmsd,rmsf,rg,energy,interactions,logs,metadata}
LOG="$OUT/logs/analysis.log"
exec > >(tee -a "$LOG") 2>&1

die() { echo "ERROR: $*" >&2; exit 1; }
need_file() { [[ -s "$1" ]] || die "Required file not found or empty: $1"; }
command -v "$GMX" >/dev/null 2>&1 || die "Cannot find GROMACS command '$GMX'. Set GMX=gmx or GMX=gmx_mpi."
need_file "$REF_TPR"
need_file "$TRAJ"
need_file "$ENERGY"
need_file "$STRUCTURE"

echo "=== GROMACS protein–ligand analysis ==="
date
echo "GROMACS command: $GMX"
echo "Reference TPR: $REF_TPR"
echo "Trajectory: $TRAJ"
echo "Energy: $ENERGY"
echo "Structure: $STRUCTURE"
echo "Ligand residue name: $LIGAND_RESNAME"
echo

# Save reproducibility metadata and validate the input files.
"$GMX" --version > "$OUT/metadata/gromacs_version.txt" 2>&1 || true
cp -f "$REF_TPR" "$OUT/metadata/" 2>/dev/null || true
"$GMX" check -f "$TRAJ" > "$OUT/logs/trajectory_check.txt" 2>&1 || {
    cat "$OUT/logs/trajectory_check.txt"
    die "GROMACS could not validate trajectory $TRAJ."
}
"$GMX" check -e "$ENERGY" > "$OUT/logs/energy_check.txt" 2>&1 || true

# Use the supplied final trajectory directly. Do not concatenate all checkpoint
# XTCs automatically: they may be cumulative and concatenating them could duplicate frames.
cp -f "$TRAJ" "$OUT/trajectories/trajectory_input.xtc"

# Build the ligand group first, then rename it and create the complex group using
# the actual group number. This avoids assuming a fixed number of default groups.
printf 'r %s\nq\n' "$LIGAND_RESNAME" \
    | "$GMX" make_ndx -f "$STRUCTURE" -o "$OUT/index/analysis_stage1.ndx" \
      > "$OUT/logs/make_ndx_stage1.log" 2>&1 || {
        cat "$OUT/logs/make_ndx_stage1.log"
        die "Could not create a ligand index group. Check LIGAND_RESNAME and the structure's residue names."
      }

python3 - "$OUT/index/analysis_stage1.ndx" "$OUT/metadata/ligand_group_id.txt" <<'PY'
import re, sys
idx, out = sys.argv[1:]
groups = []
with open(idx) as f:
    for line in f:
        m = re.match(r'\s*\[\s*(.*?)\s*\]\s*$', line)
        if m: groups.append(m.group(1))
# make_ndx appends the newly selected residue group; locate the final matching name.
hits = [i for i, name in enumerate(groups) if name not in ("System", "Protein", "Protein-H", "C-alpha", "Backbone", "MainChain", "MainChain+H", "SideChain", "SideChain-H", "Prot-Masses", "non-Protein", "Water", "SOL")]
if not hits:
    raise SystemExit(f"No custom ligand group found in {idx}. Groups: {groups}")
with open(out, "w") as f: f.write(str(hits[-1]))
print("Index groups:", list(enumerate(groups)))
PY
LIGAND_RAW_GROUP="$(cat "$OUT/metadata/ligand_group_id.txt")"
printf 'name %s Ligand\n1 | %s\nq\n' "$LIGAND_RAW_GROUP" "$LIGAND_RAW_GROUP" \
    | "$GMX" make_ndx -f "$STRUCTURE" -n "$OUT/index/analysis_stage1.ndx" \
      -o "$OUT/index/analysis.ndx" > "$OUT/logs/make_ndx.log" 2>&1 || {
        cat "$OUT/logs/make_ndx.log"
        die "Could not finalize the index groups."
      }

# Locate the final named groups from the index file rather than assuming IDs.
python3 - "$OUT/index/analysis.ndx" "$OUT/metadata/group_ids.env" <<'PY'
import re, sys
idx, out = sys.argv[1:]
groups = []
with open(idx) as f:
    for line in f:
        m = re.match(r'\s*\[\s*(.*?)\s*\]\s*$', line)
        if m: groups.append(m.group(1))
def find(name):
    hits = [i for i, g in enumerate(groups) if g == name]
    if not hits: raise SystemExit(f"Required index group {name!r} not found. Groups: {groups}")
    return hits[-1]
protein, ligand, complex_, backbone = find("Protein"), find("Ligand"), find("Protein_Ligand"), find("Backbone")
with open(out, "w") as f:
    f.write(f"PROTEIN_GROUP={protein}\nLIGAND_GROUP={ligand}\nCOMPLEX_GROUP={complex_}\nBACKBONE_GROUP={backbone}\n")
print("Index groups:")
for i, g in enumerate(groups): print(f"  {i}: {g}")
PY
# shellcheck disable=SC1090
source "$OUT/metadata/group_ids.env"

echo
echo "Selected group IDs: Protein=$PROTEIN_GROUP, Ligand=$LIGAND_GROUP, Complex=$COMPLEX_GROUP, Backbone=$BACKBONE_GROUP"
echo "Review $OUT/logs/make_ndx.log and $OUT/index/analysis.ndx if these are unexpected."

# Create a centered, whole-molecule trajectory. Center on the protein and output
# the complete system so ligand, protein, water and ions remain available.
# Interactive selections are generated from detected group IDs.
printf '%s\n%s\n' "$PROTEIN_GROUP" 0 \
  | "$GMX" trjconv -s "$REF_TPR" -f "$OUT/trajectories/trajectory_input.xtc" \
      -n "$OUT/index/analysis.ndx" -o "$OUT/trajectories/trajectory_centered.xtc" \
      -pbc mol -center -ur compact \
      > "$OUT/logs/trjconv_center.log" 2>&1 || {
        cat "$OUT/logs/trjconv_center.log"
        die "Trajectory centering failed."
      }

# Align protein backbone to the reference structure. This trajectory is intended
# for protein RMSD/RMSF/Rg and subsequent protein-relative analyses.
printf '%s\n%s\n' "$BACKBONE_GROUP" 0 \
  | "$GMX" trjconv -s "$REF_TPR" -f "$OUT/trajectories/trajectory_centered.xtc" \
      -n "$OUT/index/analysis.ndx" -o "$OUT/trajectories/trajectory_fit.xtc" \
      -fit rot+trans \
      > "$OUT/logs/trjconv_fit.log" 2>&1 || {
        cat "$OUT/logs/trjconv_fit.log"
        die "Trajectory fitting failed."
      }

# Save a final-frame coordinate file for visualization.
printf '%s\n' 0 \
  | "$GMX" trjconv -s "$REF_TPR" -f "$OUT/trajectories/trajectory_fit.xtc" \
      -n "$OUT/index/analysis.ndx" -o "$OUT/structures/final_frame.gro" -dump 500000 \
      > "$OUT/logs/final_frame.log" 2>&1 || {
        echo "WARNING: Could not extract exactly 500 ns; extracting the final frame instead."
        printf '%s\n' 0 | "$GMX" trjconv -s "$REF_TPR" \
          -f "$OUT/trajectories/trajectory_fit.xtc" -n "$OUT/index/analysis.ndx" \
          -o "$OUT/structures/final_frame.gro" \
          > "$OUT/logs/final_frame_fallback.log" 2>&1 || true
      }

# Protein backbone RMSD: fit and calculate on the backbone.
printf '%s\n%s\n' "$BACKBONE_GROUP" "$BACKBONE_GROUP" \
  | "$GMX" rms -s "$REF_TPR" -f "$OUT/trajectories/trajectory_centered.xtc" \
      -n "$OUT/index/analysis.ndx" -o "$OUT/rmsd/protein_backbone_rmsd.xvg" -tu ns \
      > "$OUT/logs/rmsd.log" 2>&1

# Ligand RMSD relative to the protein-aligned reference. The input trajectory
# is fitted to the protein backbone; select Backbone as the fit group and Ligand
# as the RMSD group when prompted.
printf '%s\n%s\n' "$BACKBONE_GROUP" "$LIGAND_GROUP" \
  | "$GMX" rms -s "$REF_TPR" -f "$OUT/trajectories/trajectory_fit.xtc" \
      -n "$OUT/index/analysis.ndx" -o "$OUT/rmsd/ligand_rmsd_protein_aligned.xvg" -tu ns \
      > "$OUT/logs/ligand_rmsd.log" 2>&1 || {
        echo "WARNING: Ligand RMSD failed; inspect $OUT/logs/ligand_rmsd.log."
      }

# Protein residue RMSF.
printf '%s\n' "$PROTEIN_GROUP" \
  | "$GMX" rmsf -s "$REF_TPR" -f "$OUT/trajectories/trajectory_fit.xtc" \
      -n "$OUT/index/analysis.ndx" -o "$OUT/rmsf/protein_rmsf.xvg" -res \
      > "$OUT/logs/rmsf.log" 2>&1

# Protein radius of gyration.
printf '%s\n' "$PROTEIN_GROUP" \
  | "$GMX" gyrate -s "$REF_TPR" -f "$OUT/trajectories/trajectory_fit.xtc" \
      -n "$OUT/index/analysis.ndx" -o "$OUT/rg/protein_rg.xvg" \
      > "$OUT/logs/gyrate.log" 2>&1

# Protein-ligand hydrogen bonds. The two group selections are prompted in order.
printf '%s\n%s\n' "$PROTEIN_GROUP" "$LIGAND_GROUP" \
  | "$GMX" hbond -s "$REF_TPR" -f "$OUT/trajectories/trajectory_fit.xtc" \
      -n "$OUT/index/analysis.ndx" -num "$OUT/interactions/hbond_count.xvg" \
      -hbn "$OUT/interactions/hbond_pairs.ndx" \
      > "$OUT/logs/hbond.log" 2>&1 || {
        echo "WARNING: gmx hbond failed; inspect $OUT/logs/hbond.log."
      }

# Protein-ligand minimum distance over time.
printf '%s\n%s\n' "$PROTEIN_GROUP" "$LIGAND_GROUP" \
  | "$GMX" mindist -s "$REF_TPR" -f "$OUT/trajectories/trajectory_fit.xtc" \
      -n "$OUT/index/analysis.ndx" -od "$OUT/interactions/protein_ligand_min_distance.xvg" \
      -on "$OUT/interactions/protein_ligand_contact_count.xvg" -d "$CUTOFF_NM" \
      > "$OUT/logs/mindist.log" 2>&1 || {
        echo "WARNING: gmx mindist failed; inspect $OUT/logs/mindist.log."
      }

# Energy extraction: interactively select common terms that exist in the EDR.
# GROMACS versions/MDP settings can change the available list, so this stage is
# interactive by design rather than guessing energy-term numbers.
echo
echo "=== Energy extraction ==="
echo "GROMACS will list energy terms from $ENERGY."
echo "At the prompt, select available terms (e.g. Potential, Kinetic, Total Energy,"
echo "Temperature, Pressure, Density) by their displayed numbers, then enter 0."
echo "Output: $OUT/energy/energy_timeseries.xvg"
"$GMX" energy -f "$ENERGY" -o "$OUT/energy/energy_timeseries.xvg" \
  2>&1 | tee "$OUT/logs/energy_extract.log" || {
    echo "WARNING: Energy extraction was not completed; rerun gmx energy interactively."
  }

# Convert XVG files to CSV and create standard plots if Python packages are available.
if python3 - <<'PY' >/dev/null 2>&1
import numpy, pandas, matplotlib
PY
then
  python3 - "$OUT" <<'PY'
import os, re, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = sys.argv[1]
xvg_files = [
    ("rmsd/protein_backbone_rmsd.xvg", "Protein backbone RMSD", "RMSD (nm)", "rmsd/protein_backbone_rmsd"),
    ("rmsf/protein_rmsf.xvg", "Protein RMSF", "RMSF (nm)", "rmsf/protein_rmsf"),
    ("rmsd/ligand_rmsd_protein_aligned.xvg", "Ligand RMSD (protein-aligned)", "RMSD (nm)", "rmsd/ligand_rmsd_protein_aligned"),
    ("rg/protein_rg.xvg", "Protein radius of gyration", "Rg (nm)", "rg/protein_rg"),
    ("interactions/hbond_count.xvg", "Protein–ligand hydrogen bonds", "Hydrogen bonds", "interactions/hbond_count"),
    ("interactions/protein_ligand_min_distance.xvg", "Protein–ligand minimum distance", "Distance (nm)", "interactions/protein_ligand_min_distance"),
    ("interactions/protein_ligand_contact_count.xvg", "Protein–ligand contacts", "Contact count", "interactions/protein_ligand_contact_count"),
    ("energy/energy_timeseries.xvg", "Energy/thermodynamic time series", "Value (term-dependent)", "energy/energy_timeseries"),
]
def read_xvg(path):
    legends, rows = [], []
    with open(path, errors="replace") as f:
        for line in f:
            if line.startswith("@"):
                m = re.search(r'@ s\d+ legend "([^"]+)"', line)
                if m: legends.append(m.group(1))
            elif line.strip() and not line.startswith(("#", "&")):
                try: rows.append([float(x) for x in line.split()])
                except ValueError: pass
    if not rows: return None, legends
    n = min(map(len, rows))
    return np.asarray([r[:n] for r in rows]), legends

for rel, title, ylabel, stem in xvg_files:
    path = os.path.join(root, rel)
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        continue
    arr, legends = read_xvg(path)
    if arr is None: continue
    cols = ["time"] + [legends[i] if i < len(legends) else f"series_{i+1}" for i in range(arr.shape[1]-1)]
    # RMSF x-axis is residue index, not time.
    if rel.startswith("rmsf/"): cols[0] = "residue_index"
    df = pd.DataFrame(arr, columns=cols)
    csv_path = os.path.join(root, stem + ".csv")
    df.to_csv(csv_path, index=False)
    fig, ax = plt.subplots(figsize=(9, 5))
    x = arr[:, 0]
    for j in range(1, arr.shape[1]):
        label = cols[j]
        ax.plot(x, arr[:, j], label=label)
    ax.set_title(title)
    ax.set_xlabel("Residue index" if rel.startswith("rmsf/") else "Time (ns)")
    ax.set_ylabel(ylabel)
    if arr.shape[1] > 2: ax.legend(fontsize="small")
    fig.tight_layout()
    fig.savefig(os.path.join(root, stem + ".png"), dpi=180)
    plt.close(fig)
print("Converted available XVG files to CSV and PNG.")
PY
else
  echo "Python plotting dependencies not available; XVG outputs were retained. Install numpy pandas matplotlib to create CSV/PNG plots."
fi

# Summary of input files and output inventory.
{
  echo "Analysis date: $(date -Is)"
  echo "GROMACS: $GMX"
  echo "Reference TPR: $REF_TPR"
  echo "Input trajectory: $TRAJ"
  echo "Energy file: $ENERGY"
  echo "Structure used for index: $STRUCTURE"
  echo "Ligand residue name: $LIGAND_RESNAME"
  echo "Contact cutoff (nm): $CUTOFF_NM"
  echo
  echo "Trajectory check:"
  cat "$OUT/logs/trajectory_check.txt"
} > "$OUT/metadata/summary.txt"

echo
echo "=== Analysis pipeline finished ==="
echo "Results are in: $OUT/"
echo "Read first: $OUT/metadata/summary.txt"
echo "Check logs in: $OUT/logs/"
echo "NOTE: per-residue contact occupancy is not computed by this GROMACS-only script; contact count and minimum distance are included."
echo "For those, use a protein-aligned ligand trajectory and an MDAnalysis post-processing script."

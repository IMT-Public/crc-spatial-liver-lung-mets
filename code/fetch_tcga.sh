#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Fetch the public TCGA-COAD/READ files that 04_figure4.R needs for Fig 4d.
#
#   cd code && ./fetch_tcga.sh
#
# Source: cBioPortal study `coadread_tcga_pan_can_atlas_2018`, served from the
# cBioPortal datahub repository. The two files are stored there with Git LFS,
# so they are fetched through the media.githubusercontent.com LFS endpoint
# rather than raw.githubusercontent.com, which would return the LFS pointer.
#
# To download by hand instead, use the study page:
#   https://www.cbioportal.org/study/summary?id=coadread_tcga_pan_can_atlas_2018
#
# These files are public and redistributable, so they are BUNDLED into /data
# for the Code Ocean capsule: a Reproducible Run has no network access. This
# script is for local and CI runs. Neither file is tracked in git -- the
# expression matrix alone is 83 MB.
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "$0")"

DEST="${CRC_DATA_DIR:-../data}/data_tcga"
STUDY="coadread_tcga_pan_can_atlas_2018"
BASE="https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/${STUDY}"

FILES=(
  data_clinical_patient.txt
  data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt
)

mkdir -p "$DEST"

for f in "${FILES[@]}"; do
  if [ -s "$DEST/$f" ]; then
    echo "already present: $DEST/$f"
    continue
  fi
  echo "Downloading $f from the cBioPortal datahub..."
  tmp="$DEST/.$f.part"
  curl -fL --progress-bar "$BASE/$f" -o "$tmp"

  # A Git LFS pointer is ~130 bytes and starts with "version https://git-lfs".
  # Catching that here beats failing three minutes into 04_figure4.R.
  if head -c 24 "$tmp" | grep -q '^version https://git-lfs'; then
    rm -f "$tmp"
    echo "ERROR: got an LFS pointer, not the file. The datahub layout may have" >&2
    echo "       changed -- download from the study page instead:" >&2
    echo "       https://www.cbioportal.org/study/summary?id=${STUDY}" >&2
    exit 1
  fi

  mv "$tmp" "$DEST/$f"
  echo "wrote $DEST/$f ($(wc -c < "$DEST/$f") bytes)"
done

echo
echo "Done. The third file 04_figure4.R needs, some_cell_subtype_markers.xlsx,"
echo "is specific to this study and ships with the data -- see data/README.md."

#!/usr/bin/env bash
# Controlled Fishbone--Moncrief Kerr--Schild spin sweep.
#
# The base setup is preserved; only a_spin and output_folder are changed.
# initialization.f90 is restored from an exact backup at the end of the sweep.

set -euo pipefail

source_file="initialization.f90"
build_jobs="${BUILD_JOBS:-4}"
omp_threads="${OMP_NUM_THREADS:-4}"
campaign_log="${CAMPAIGN_LOG:-log_FM_KS_metricfix_campaign.txt}"
python_command="${PYTHON_COMMAND:-python}"
plot_script="${PLOT_SCRIPT:-plot_fm_profiles.py}"
spins=(0.0d0 0.2d0 0.5d0 0.9d0)

if [[ ! -f "$source_file" ]]; then
  echo "ERROR: run this script from the GRHD source directory." >&2
  exit 1
fi
if [[ ! -f "$plot_script" ]]; then
  echo "ERROR: radial-profile script not found: $plot_script" >&2
  exit 1
fi
if ! MPLBACKEND=Agg "$python_command" -c 'import matplotlib, numpy, pyvista'; then
  echo "ERROR: Python plotting dependencies are not available." >&2
  exit 1
fi

backup_file="$(mktemp ./initialization.f90.sweep.XXXXXX)"
diff_file="$(mktemp ./initialization.f90.diff.XXXXXX)"
cp -p -- "$source_file" "$backup_file"

restore_source() {
  if [[ -e "$backup_file" ]]; then
    mv -f -- "$backup_file" "$source_file"
  fi
  rm -f -- "$diff_file"
}

on_interrupt() {
  local signal_number="$1"
  exit "$((128 + signal_number))"
}

validate_temporary_diff() {
  local status=0

  diff -u --label initialization.f90.before --label initialization.f90.temporary \
    "$backup_file" "$source_file" > "$diff_file" || status=$?
  if (( status > 1 )); then
    echo "ERROR: could not inspect the temporary initialization.f90 diff." >&2
    return 1
  fi

  if ! awk '
    /^--- / || /^\+\+\+ / || /^@@/ || /^[^+-]/ { next }
    {
      line = substr($0, 2)
      if (line !~ /^[[:space:]]*a_spin[[:space:]]*=/ &&
          line !~ /^[[:space:]]*output_folder[[:space:]]*=/) {
        print "UNAUTHORIZED TEMPORARY DIFF: " $0 > "/dev/stderr"
        bad = 1
      }
    }
    END { exit bad }
  ' "$diff_file"; then
    echo "ERROR: initialization.f90 changed outside a_spin/output_folder; aborting." >&2
    return 1
  fi
}

trap restore_source EXIT
trap 'on_interrupt 2' INT
trap 'on_interrupt 15' TERM

if [[ -e "$campaign_log" ]]; then
  echo "ERROR: campaign log already exists: $campaign_log" >&2
  exit 1
fi

for spin in "${spins[@]}"; do
  spin_tag="${spin%.0d0}"
  spin_tag="${spin_tag%d0}"
  prefix="FM_KS_metricfix_a${spin_tag}"
  folder="${prefix}_data"
  log_file="log_${prefix}.txt"
  plot_file="${prefix}_perfil_radial.png"

  if [[ -e "$folder" || -e "$log_file" || -e "$plot_file" ]]; then
    echo "ERROR: refusing to overwrite existing output, log, or plot for ${prefix}." >&2
    exit 1
  fi
done

exec > >(tee -- "$campaign_log") 2>&1

for spin in "${spins[@]}"; do
  spin_tag="${spin%.0d0}"
  spin_tag="${spin_tag%d0}"
  prefix="FM_KS_metricfix_a${spin_tag}"
  folder="${prefix}_data"
  log_file="log_${prefix}.txt"
  plot_file="${prefix}_perfil_radial.png"

  cp -p -- "$backup_file" "$source_file"
  sed -i -E \
    -e "/subroutine setup_fishbone_moncrief_equatorial\(\)/,/end subroutine setup_fishbone_moncrief_equatorial/ s/^[[:space:]]*a_spin[[:space:]]*=.*/    a_spin = ${spin}/" \
    -e "/subroutine setup_fishbone_moncrief_equatorial\(\)/,/end subroutine setup_fishbone_moncrief_equatorial/ s/^[[:space:]]*output_folder[[:space:]]*=.*/    output_folder = '${folder}'/" \
    "$source_file"

  validate_temporary_diff

  if [[ "$(sed -n '/subroutine setup_fishbone_moncrief_equatorial()/,/end subroutine setup_fishbone_moncrief_equatorial/ { /^[[:space:]]*a_spin[[:space:]]*=/p; }' "$source_file")" != "    a_spin = ${spin}" ]]; then
    echo "ERROR: a_spin was not set exactly once inside the FM equatorial setup." >&2
    exit 1
  fi
  if [[ "$(sed -n '/subroutine setup_fishbone_moncrief_equatorial()/,/end subroutine setup_fishbone_moncrief_equatorial/ { /^[[:space:]]*output_folder[[:space:]]*=/p; }' "$source_file")" != "    output_folder = '${folder}'" ]]; then
    echo "ERROR: output_folder was not set exactly once inside the FM equatorial setup." >&2
    exit 1
  fi

  echo "=== ${prefix}: starting with OMP_NUM_THREADS=${omp_threads} ==="
  make clean
  make -j"${build_jobs}"
  OMP_NUM_THREADS="$omp_threads" ./grhd2 2>&1 | tee -- "$log_file"

  MPLBACKEND=Agg "$python_command" "$plot_script" \
    "${folder}/*/*.vtk" \
    --nx 400 --frames 1000 --xlim 1 40 --ylim -0.05 1.2 \
    --output "$plot_file" --no-show
  if [[ ! -s "$plot_file" ]]; then
    echo "ERROR: radial-profile plot was not created: $plot_file" >&2
    exit 1
  fi

  echo "=== ${prefix}: radial profile saved to ${plot_file} ==="
  echo "=== ${prefix} completed ==="
done

echo "All corrected Kerr--Schild spin runs completed."

#!/usr/bin/env bash
# Controlled Fishbone--Moncrief Kerr--Schild spin sweep.
#
# It preserves the source tree: initialization.f90 is restored on exit even
# if the session is interrupted.  r_min=2 reproduces the legacy campaign so
# any difference can be attributed to the corrected metric/inversion code.

set -euo pipefail

source_file="initialization.f90"
build_jobs="${BUILD_JOBS:-4}"
omp_threads="${OMP_NUM_THREADS:-4}"
spins=(0.0d0 0.2d0 0.5d0 0.9d0)
legacy_r_min="2.0d0"

if [[ ! -f "$source_file" ]]; then
  echo "ERROR: run this script from the GRHD source directory." >&2
  exit 1
fi

backup_file="$(mktemp ./initialization.f90.sweep.XXXXXX)"
cp "$source_file" "$backup_file"

restore_source() {
  cp "$backup_file" "$source_file"
  rm -f "$backup_file"
}
trap restore_source EXIT INT TERM

for spin in "${spins[@]}"; do
  spin_tag="${spin%.0d0}"
  spin_tag="${spin_tag%d0}"
  prefix="FM_KS_metricfix_a${spin_tag}"
  folder="${prefix}_data"
  log_file="log_${prefix}.txt"

  cp "$backup_file" "$source_file"
  sed -i -E \
    -e "/subroutine setup_fishbone_moncrief_equatorial\(\)/,/end subroutine setup_fishbone_moncrief_equatorial/ s/^[[:space:]]*a_spin[[:space:]]*=.*/    a_spin = ${spin}/" \
    -e "/subroutine setup_fishbone_moncrief_equatorial\(\)/,/end subroutine setup_fishbone_moncrief_equatorial/ s/^[[:space:]]*nx[[:space:]]*=.*/    nx = 400 ; r_min = ${legacy_r_min} ; r_max = 40.0d0/" \
    -e "/subroutine setup_fishbone_moncrief_equatorial\(\)/,/end subroutine setup_fishbone_moncrief_equatorial/ s/^[[:space:]]*final_time[[:space:]]*=.*/    final_time = 1000.0d0/" \
    -e "/subroutine setup_fishbone_moncrief_equatorial\(\)/,/end subroutine setup_fishbone_moncrief_equatorial/ s/^[[:space:]]*output_prefix[[:space:]]*=.*/    output_prefix = '${prefix}'/" \
    -e "/subroutine setup_fishbone_moncrief_equatorial\(\)/,/end subroutine setup_fishbone_moncrief_equatorial/ s/^[[:space:]]*output_folder[[:space:]]*=.*/    output_folder = '${folder}'/" \
    "$source_file"

  echo "=== ${prefix}: T=1000, r_min=${legacy_r_min}, OMP=${omp_threads} ==="
  make clean
  make -j"${build_jobs}"
  OMP_NUM_THREADS="$omp_threads" ./grhd2 > "$log_file" 2>&1
  echo "=== ${prefix} completed ==="
done

echo "All corrected Kerr--Schild spin runs completed."

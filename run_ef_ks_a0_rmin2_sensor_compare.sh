#!/usr/bin/env bash
# Controlled Fishbone--Moncrief comparison between EF and KS at a=0.
#
# Both runs use the same logarithmic grid, r_min=2, shock sensor, WENO5/HLLE,
# and T=1000. The base initialization file is restored exactly on every exit.

set -euo pipefail

source_file="initialization.f90"
build_jobs="${BUILD_JOBS:-4}"
omp_threads="${OMP_NUM_THREADS:-4}"
python_command="${PYTHON_COMMAND:-python}"
plot_script="${PLOT_SCRIPT:-plot_fm_profiles.py}"
campaign_log="${CAMPAIGN_LOG:-log_FM_EFKS_compare_a0_rmin2_sensor_campaign.txt}"

metric_tags=(EF KS)
metric_names=(Eddington-Finkelstein Kerr-Schild)

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

backup_file="$(mktemp ./initialization.f90.compare.XXXXXX)"
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
      sub(/\r$/, "", line)
      if (line ~ /^[[:space:]]*metric_type[[:space:]]*=/ ||
          line ~ /^[[:space:]]*use_shock_sensor[[:space:]]*=/ ||
          line ~ /^[[:space:]]*a_spin[[:space:]]*=/ ||
          line ~ /^[[:space:]]*output_prefix[[:space:]]*=/ ||
          line ~ /^[[:space:]]*output_folder[[:space:]]*=/ ||
          line ~ /^[[:space:]]*nx[[:space:]]*=[[:space:]]*400[[:space:]]*;[[:space:]]*r_min[[:space:]]*=[[:space:]]*(1\.2d0|2\.0d0)[[:space:]]*;[[:space:]]*r_max[[:space:]]*=[[:space:]]*40\.0d0[[:space:]]*$/) {
        next
      }
      print "UNAUTHORIZED TEMPORARY DIFF: " $0 > "/dev/stderr"
      bad = 1
    }
    END { exit bad }
  ' "$diff_file"; then
    echo "ERROR: initialization.f90 changed outside the comparison parameters; aborting." >&2
    return 1
  fi
}

extract_fm_assignment() {
  local variable="$1"
  sed -n "/subroutine setup_fishbone_moncrief_equatorial()/,/end subroutine setup_fishbone_moncrief_equatorial/ { /^[[:space:]]*${variable}[[:space:]]*=/p; }" "$source_file" | tr -d '\r'
}

trap restore_source EXIT
trap 'on_interrupt 1' HUP
trap 'on_interrupt 2' INT
trap 'on_interrupt 15' TERM

if [[ -e "$campaign_log" ]]; then
  echo "ERROR: campaign log already exists: $campaign_log" >&2
  exit 1
fi

for tag in "${metric_tags[@]}"; do
  prefix="FM_EFKS_compare_${tag}_a0_rmin2_sensor"
  folder="${prefix}_data"
  log_file="log_${prefix}.txt"
  plot_file="perfil_radial_${prefix}_T1000.png"

  if [[ -e "$folder" || -e "$log_file" || -e "$plot_file" ]]; then
    echo "ERROR: refusing to overwrite existing output, log, or plot for ${prefix}." >&2
    exit 1
  fi
done

exec > >(tee -- "$campaign_log") 2>&1

for index in "${!metric_tags[@]}"; do
  tag="${metric_tags[$index]}"
  metric="${metric_names[$index]}"
  prefix="FM_EFKS_compare_${tag}_a0_rmin2_sensor"
  folder="${prefix}_data"
  log_file="log_${prefix}.txt"
  plot_file="perfil_radial_${prefix}_T1000.png"

  cp -p -- "$backup_file" "$source_file"
  sed -i -E \
    -e "/subroutine setup_fishbone_moncrief_equatorial\(\)/,/end subroutine setup_fishbone_moncrief_equatorial/ s/^[[:space:]]*metric_type[[:space:]]*=.*/    metric_type = '${metric}'/" \
    -e "/subroutine setup_fishbone_moncrief_equatorial\(\)/,/end subroutine setup_fishbone_moncrief_equatorial/ s/^[[:space:]]*use_shock_sensor[[:space:]]*=.*/    use_shock_sensor = .true./" \
    -e "/subroutine setup_fishbone_moncrief_equatorial\(\)/,/end subroutine setup_fishbone_moncrief_equatorial/ s/^[[:space:]]*a_spin[[:space:]]*=.*/    a_spin = 0.0d0/" \
    -e "/subroutine setup_fishbone_moncrief_equatorial\(\)/,/end subroutine setup_fishbone_moncrief_equatorial/ s/(r_min[[:space:]]*=[[:space:]]*)[^;]*/\12.0d0 /" \
    -e "/subroutine setup_fishbone_moncrief_equatorial\(\)/,/end subroutine setup_fishbone_moncrief_equatorial/ s/^[[:space:]]*output_prefix[[:space:]]*=.*/    output_prefix = '${prefix}'/" \
    -e "/subroutine setup_fishbone_moncrief_equatorial\(\)/,/end subroutine setup_fishbone_moncrief_equatorial/ s/^[[:space:]]*output_folder[[:space:]]*=.*/    output_folder = '${folder}'/" \
    "$source_file"

  validate_temporary_diff

  [[ "$(extract_fm_assignment metric_type)" == "    metric_type = '${metric}'" ]]
  [[ "$(extract_fm_assignment use_shock_sensor)" == "    use_shock_sensor = .true." ]]
  [[ "$(extract_fm_assignment a_spin)" == "    a_spin = 0.0d0" ]]
  [[ "$(extract_fm_assignment nx)" == "    nx = 400 ; r_min = 2.0d0 ; r_max = 40.0d0" ]]
  [[ "$(extract_fm_assignment output_prefix)" == "    output_prefix = '${prefix}'" ]]
  [[ "$(extract_fm_assignment output_folder)" == "    output_folder = '${folder}'" ]]

  echo "=== ${prefix}: starting with OMP_NUM_THREADS=${omp_threads} ==="
  make clean
  make -j"${build_jobs}"
  OMP_NUM_THREADS="$omp_threads" ./grhd2 2>&1 | tee -- "$log_file"

  MPLBACKEND=Agg "$python_command" "$plot_script" \
    "${folder}/*/*.vtk" \
    --nx 400 --frames 1000 --xlim 2 40 --ylim -0.05 1.2 \
    --output "$plot_file" --no-show
  if [[ ! -s "$plot_file" ]]; then
    echo "ERROR: radial-profile plot was not created: $plot_file" >&2
    exit 1
  fi

  echo "=== ${prefix}: radial profile saved to ${plot_file} ==="
  echo "=== ${prefix} completed ==="
done

echo "EF/KS a=0 r_min=2 shock-sensor comparison completed."

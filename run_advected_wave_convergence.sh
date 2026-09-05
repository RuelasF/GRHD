#!/usr/bin/env bash
set -euo pipefail

source_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
result_name=${1:-AdvectedWave_HLLC_Convergence_$(date +%Y%m%d_%H%M%S)}

if [[ ${result_name} = /* ]]; then
  result_root=${result_name}
else
  result_root=${source_dir}/${result_name}
fi

if [[ -e ${result_root} ]]; then
  printf 'ERROR: result directory already exists: %s\n' "${result_root}" >&2
  exit 1
fi

mkdir -p "${result_root}/data" "${result_root}/logs" \
  "${result_root}/.matplotlib-cache"

make -C "${source_dir}" -B

resolutions=(32 64 128 256 512 1024)
methods=(
  "1:godunov"
  "2:mc"
  "3:weno3"
  "4:mp5"
  "5:weno5"
)

run_method() {
  local specification=$1
  local reconstruction_id=${specification%%:*}
  local method_name=${specification#*:}
  local resolution
  local log_path

  for resolution in "${resolutions[@]}"; do
    log_path=${result_root}/logs/${method_name}_N${resolution}.log
    (
      cd "${result_root}"
      env \
        OMP_NUM_THREADS=1 \
        GRHD_CASE_ID=6 \
        GRHD_RECONSTRUCTION="${reconstruction_id}" \
        GRHD_TVD_LIMITER=3 \
        GRHD_RIEMANN_SOLVER=2 \
        GRHD_NX="${resolution}" \
        GRHD_FINAL_TIME=0.7 \
        GRHD_CFL=0.1 \
        GRHD_SAVE_INTERVAL=10.0 \
        GRHD_OUTPUT_FOLDER=data \
        "${source_dir}/grhd2" > "${log_path}" 2>&1
    )
    if ! grep -q 'SIMULACIÓN FINALIZADA EXITOSAMENTE' "${log_path}"; then
      printf 'ERROR: failed run %s N=%s; see %s\n' "${method_name}" "${resolution}" "${log_path}" >&2
      return 1
    fi
    printf 'Completed %-8s N=%s\n' "${method_name}" "${resolution}"
  done
}

printf '%s\n' \
  'case=advected_wave' \
  'riemann_solver=HLLC' \
  'resolutions=32,64,128,256,512,1024' \
  'reconstructions=Godunov,TVD-MC,WENO3,MP5,WENO5' \
  'final_time=0.7' \
  'CFL=0.1' \
  'time_integrator=SSP-RK3' \
  'shock_sensor=false' > "${result_root}/campaign_manifest.txt"

if [[ ${GRHD_CONVERGENCE_PARALLEL:-0} == 1 ]]; then
  pids=()
  for specification in "${methods[@]}"; do
    run_method "${specification}" &
    pids+=("$!")
  done
  for pid in "${pids[@]}"; do
    wait "${pid}"
  done
else
  for specification in "${methods[@]}"; do
    run_method "${specification}"
  done
fi

env MPLCONFIGDIR="${result_root}/.matplotlib-cache" \
  python3 "${source_dir}/plot_advected_wave_convergence.py" \
  "${result_root}/data" --output-directory "${result_root}" \
  > "${result_root}/convergence_summary.txt"

printf 'Results: %s\n' "${result_root}"
printf 'Norm plot: %s\n' "${result_root}/advected_wave_hllc_norms.png"
printf 'Order plot: %s\n' "${result_root}/advected_wave_hllc_orders.png"

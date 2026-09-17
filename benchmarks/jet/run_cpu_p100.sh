#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  printf '%s\n' \
    'Uso: run_cpu_p100.sh [--cpu-only | --gpu-only]' \
    '' \
    'Variables opcionales:' \
    '  CPU_THREADS=N       Hilos OpenMP para la referencia CPU (default: CPUs visibles).' \
    '  GPU_ID=N            GPU CUDA utilizada (default: 0).' \
    '  BENCHMARK_DIR=RUTA  Directorio de resultados.' \
    '  REL_TOL=VALOR       Tolerancia L2 relativa (default: 1e-8).' \
    '  ABS_TOL=VALOR       Tolerancia Linf absoluta (default: 1e-12).'
}

run_cpu=true
run_gpu=true
case "${1:-}" in
  '') ;;
  --cpu-only) run_gpu=false ;;
  --gpu-only) run_cpu=false ;;
  -h|--help) usage; exit 0 ;;
  *) usage >&2; exit 2 ;;
esac
if (( $# > 1 )); then
  usage >&2
  exit 2
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$script_dir/../.." && pwd)
parameter_file="$script_dir/jet_p100_benchmark.par"
timestamp=$(date +%Y%m%d_%H%M%S)
result_root=${BENCHMARK_DIR:-"$repo_root/.benchmark/jet_cpu_p100_$timestamp"}
cpu_threads=${CPU_THREADS:-$(getconf _NPROCESSORS_ONLN)}
gpu_id=${GPU_ID:-0}
relative_tolerance=${REL_TOL:-1e-8}
absolute_tolerance=${ABS_TOL:-1e-12}

for command in make python3; do
  if ! command -v "$command" >/dev/null 2>&1; then
    printf 'ERROR: no se encontro %s.\n' "$command" >&2
    exit 1
  fi
done
if [[ ! "$cpu_threads" =~ ^[1-9][0-9]*$ ]]; then
  printf 'ERROR: CPU_THREADS debe ser un entero positivo.\n' >&2
  exit 1
fi
if $run_cpu && ! command -v gfortran >/dev/null 2>&1; then
  printf 'ERROR: se requiere gfortran para la referencia CPU.\n' >&2
  exit 1
fi
if $run_gpu; then
  if ! command -v nvfortran >/dev/null 2>&1; then
    printf 'ERROR: se requiere nvfortran (NVIDIA HPC SDK) para la P100.\n' >&2
    exit 1
  fi
  if ! command -v nvidia-smi >/dev/null 2>&1 || ! nvidia-smi -L >/dev/null 2>&1; then
    printf 'ERROR: nvidia-smi no puede consultar la GPU; verifique el controlador NVIDIA.\n' >&2
    exit 1
  fi
fi

mkdir -p "$result_root/bin" "$result_root/cpu" "$result_root/gpu"
metadata="$result_root/system.txt"
{
  printf 'date=%s\n' "$(date --iso-8601=seconds)"
  printf 'commit=%s\n' "$(git -C "$repo_root" rev-parse HEAD)"
  printf 'cpu_threads=%s\n' "$cpu_threads"
  printf 'gpu_id=%s\n' "$gpu_id"
  printf '\n[lscpu]\n'
  lscpu
  if command -v gfortran >/dev/null 2>&1; then
    printf '\n[gfortran]\n'
    gfortran --version || printf 'no disponible\n'
  fi
  if command -v nvfortran >/dev/null 2>&1; then
    printf '\n[nvfortran]\n'
    nvfortran --version || printf 'no disponible\n'
  fi
  if command -v nvidia-smi >/dev/null 2>&1; then
    printf '\n[nvidia-smi]\n'
    nvidia-smi --query-gpu=index,name,driver_version,memory.total --format=csv,noheader || \
      printf 'controlador o GPU no disponibles\n'
  fi
} >"$metadata"

if $run_cpu; then
  printf '[1/4] Compilando referencia OpenMP CPU...\n'
  make -C "$repo_root" clean
  make -C "$repo_root" grhd2
  cp "$repo_root/grhd2" "$result_root/bin/grhd2_cpu"
  "$result_root/bin/grhd2_cpu" --check "$parameter_file" >"$result_root/cpu/check.log" 2>&1

  printf '[2/4] Ejecutando jet CPU con %s hilos...\n' "$cpu_threads"
  cpu_start=$(python3 -c 'import time; print(time.monotonic_ns())')
  (
    cd "$result_root/cpu"
    env OMP_NUM_THREADS="$cpu_threads" OMP_DYNAMIC=FALSE OMP_PROC_BIND=close \
      OMP_PLACES=cores GFORTRAN_UNBUFFERED_ALL=y \
      "$result_root/bin/grhd2_cpu" "$parameter_file" >run.log 2>&1
  )
  cpu_end=$(python3 -c 'import time; print(time.monotonic_ns())')
  awk -v start="$cpu_start" -v finish="$cpu_end" \
    'BEGIN { printf "%.9f\n", (finish-start)/1000000000 }' \
    >"$result_root/cpu/wall_seconds.txt"
fi

if $run_gpu; then
  printf '[3/4] Compilando backend OpenACC para Tesla P100 (cc60, FP64)...\n'
  make -C "$repo_root" openacc-p100
  cp "$repo_root/grhd2" "$result_root/bin/grhd2_p100"
  CUDA_VISIBLE_DEVICES="$gpu_id" "$result_root/bin/grhd2_p100" --check "$parameter_file" \
    >"$result_root/gpu/check.log" 2>&1

  printf '[4/4] Ejecutando jet en la P100 %s...\n' "$gpu_id"
  gpu_start=$(python3 -c 'import time; print(time.monotonic_ns())')
  (
    cd "$result_root/gpu"
    env CUDA_VISIBLE_DEVICES="$gpu_id" \
      "$result_root/bin/grhd2_p100" "$parameter_file" >run.log 2>&1
  )
  gpu_end=$(python3 -c 'import time; print(time.monotonic_ns())')
  awk -v start="$gpu_start" -v finish="$gpu_end" \
    'BEGIN { printf "%.9f\n", (finish-start)/1000000000 }' \
    >"$result_root/gpu/wall_seconds.txt"
fi

make -C "$repo_root" clean >/dev/null

if $run_cpu && $run_gpu; then
  cpu_vtk="$result_root/cpu/jet_output/weno5_hlle/JET_BENCH_weno5_hlle_0001.vtk"
  gpu_vtk="$result_root/gpu/jet_output/weno5_hlle/JET_BENCH_weno5_hlle_0001.vtk"
  set +e
  python3 "$script_dir/compare_vtk.py" "$cpu_vtk" "$gpu_vtk" \
    --relative-tolerance "$relative_tolerance" \
    --absolute-tolerance "$absolute_tolerance" | tee "$result_root/comparison.txt"
  comparison_status=${PIPESTATUS[0]}
  set -e

  cpu_seconds=$(<"$result_root/cpu/wall_seconds.txt")
  gpu_seconds=$(<"$result_root/gpu/wall_seconds.txt")
  speedup=$(awk -v cpu="$cpu_seconds" -v gpu="$gpu_seconds" \
    'BEGIN { if (gpu > 0) printf "%.6f", cpu/gpu; else print "inf" }')
  {
    printf 'cpu_wall_seconds=%s\n' "$cpu_seconds"
    printf 'gpu_wall_seconds=%s\n' "$gpu_seconds"
    printf 'speedup_cpu_over_gpu=%s\n' "$speedup"
    printf 'comparison_exit_status=%s\n' "$comparison_status"
  } | tee "$result_root/timing_summary.txt"

  if (( comparison_status != 0 )); then
    printf 'ADVERTENCIA: la comparacion numerica requiere revision; consulte %s.\n' \
      "$result_root/comparison.txt" >&2
    exit "$comparison_status"
  fi
fi

printf 'Benchmark terminado. Resultados: %s\n' "$result_root"

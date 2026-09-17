#!/usr/bin/env bash
#SBATCH --job-name=inundaciones
#SBATCH --cpus-per-task=72
#PBS -N inundaciones
#PBS -l select=1:ncpus=72

# Campaña PPI ecuatorial 2.5D para EF y KS.
# Uso directo persistente: bash inundaciones.sh --detach
# Uso con Slurm:           sbatch inundaciones.sh
# Uso con PBS:             qsub inundaciones.sh

set -Eeuo pipefail
umask 027

readonly SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)/$(basename "${BASH_SOURCE[0]}")"
readonly SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
readonly SCRIPT_PARENT="$(cd "$(dirname "$SCRIPT_PATH")/.." && pwd -P)"
CAMPAIGN_ROOT="${CAMPAIGN_ROOT:-$SCRIPT_PARENT/inundaciones_ppi_800_t5000}"
CONCURRENT_CASES="${CONCURRENT_CASES:-6}"
PREFLIGHT_ONLY="${PREFLIGHT_ONLY:-0}"
PRIMARY_THREADS="${PRIMARY_THREADS:-72}"
FALLBACK_THREADS="${FALLBACK_THREADS:-60}"
RECONSTRUCTIONS="${RECONSTRUCTIONS:-weno5 mp5}"
LAUNCHER_NAME="${LAUNCHER_NAME:-$(basename "$SCRIPT_PATH")}"
CAMPAIGN_LABEL="${CAMPAIGN_LABEL:-PPI MP5/WENO5}"
CPU_BINDING="${CPU_BINDING:-none}"

if [[ "$CAMPAIGN_ROOT" =~ [[:space:]] ]]; then
  printf 'ERROR: CAMPAIGN_ROOT no puede contener espacios: %s\n' "$CAMPAIGN_ROOT" >&2
  exit 2
fi

show_help() {
  cat <<EOF
$CAMPAIGN_LABEL: malla 800x1x400, t=5000 M, perturbación en t=1000 M.

  bash $LAUNCHER_NAME --detach   Lanza con nohup y devuelve el PID.
  bash $LAUNCHER_NAME           Ejecuta en primer plano.
  sbatch $LAUNCHER_NAME         Envía a Slurm con nombre inundaciones.
  qsub $LAUNCHER_NAME           Envía a PBS con nombre inundaciones.

Variables opcionales: CAMPAIGN_ROOT, GRHD_EXECUTABLE y CPU_BINDING.
PREFLIGHT_ONLY=1 valida el ejecutable, los hilos y los .par sin correr.
Al reejecutar sobre la misma CAMPAIGN_ROOT se conserva el binario inicial,
se omiten casos terminados y los incompletos continúan desde su checkpoint.
EOF
}

case "${1:-}" in
  --help|-h)
    show_help
    exit 0
    ;;
  --detach)
    mkdir -p "$CAMPAIGN_ROOT"
    nohup bash "$SCRIPT_PATH" --run >"$CAMPAIGN_ROOT/inundaciones.bootstrap.log" 2>&1 </dev/null &
    printf 'inundaciones iniciado: PID=%s\n' "$!"
    printf 'Log inicial: %s/inundaciones.bootstrap.log\n' "$CAMPAIGN_ROOT"
    exit 0
    ;;
  --run|'')
    ;;
  *)
    printf 'ERROR: opción desconocida: %s\n' "$1" >&2
    show_help >&2
    exit 2
    ;;
esac

SOURCE_EXECUTABLE="${GRHD_EXECUTABLE:-$SCRIPT_DIR/grhd2}"
readonly PARAMETER_ROOT="$CAMPAIGN_ROOT/parameters"
readonly RESULT_ROOT="$CAMPAIGN_ROOT/results"
readonly LOG_ROOT="$CAMPAIGN_ROOT/logs"
readonly STATE_ROOT="$CAMPAIGN_ROOT/state"
readonly BIN_ROOT="$CAMPAIGN_ROOT/bin"
readonly DRIVER_LOG="$LOG_ROOT/inundaciones.driver.log"
readonly SUMMARY="$RESULT_ROOT/summary.tsv"
readonly EXECUTABLE_HASH_FILE="$STATE_ROOT/inundaciones.sha256"
readonly CAMPAIGN_EXECUTABLE="$BIN_ROOT/inundaciones"

mkdir -p "$PARAMETER_ROOT" "$RESULT_ROOT" "$LOG_ROOT" "$STATE_ROOT" "$BIN_ROOT"
exec > >(tee -a "$DRIVER_LOG") 2>&1

timestamp() {
  date -Is
}

die() {
  printf '[%s] ERROR: %s\n' "$(timestamp)" "$*" >&2
  exit 1
}

terminate_children() {
  local signal_name="$1" pid
  trap - HUP INT TERM
  printf '[%s] Señal %s; deteniendo procesos hijos. Los checkpoints se conservan.\n' \
    "$(timestamp)" "$signal_name"
  while read -r pid; do
    [[ -n "$pid" ]] && kill -TERM "$pid" 2>/dev/null || true
  done < <(jobs -pr)
  wait || true
  exit 130
}
trap 'terminate_children HUP' HUP
trap 'terminate_children INT' INT
trap 'terminate_children TERM' TERM

for command_name in awk sort find sha256sum nproc; do
  command -v "$command_name" >/dev/null 2>&1 || die "falta el comando requerido: $command_name"
done

prepare_executable() {
  local source_hash campaign_hash pinned_hash

  [[ -f "$SOURCE_EXECUTABLE" ]] || die "no se encontró el ejecutable: $SOURCE_EXECUTABLE"
  [[ -x "$SOURCE_EXECUTABLE" ]] || chmod u+x "$SOURCE_EXECUTABLE" || \
    die "no se pudo marcar como ejecutable: $SOURCE_EXECUTABLE"
  source_hash="$(sha256sum "$SOURCE_EXECUTABLE" | awk '{print $1}')"

  if [[ ! -x "$CAMPAIGN_EXECUTABLE" ]]; then
    cp -p "$SOURCE_EXECUTABLE" "$CAMPAIGN_EXECUTABLE"
    chmod 750 "$CAMPAIGN_EXECUTABLE"
  fi
  campaign_hash="$(sha256sum "$CAMPAIGN_EXECUTABLE" | awk '{print $1}')"
  [[ "$campaign_hash" == "$source_hash" ]] || \
    die 'el ejecutable recibido no coincide con el binario fijado para esta campaña'

  if [[ -s "$EXECUTABLE_HASH_FILE" ]]; then
    pinned_hash="$(awk '{print $1}' "$EXECUTABLE_HASH_FILE")"
    [[ "$campaign_hash" == "$pinned_hash" ]] || \
      die 'cambió el hash del ejecutable; use otra CAMPAIGN_ROOT para una campaña nueva'
  else
    printf '%s  inundaciones\n' "$campaign_hash" >"$EXECUTABLE_HASH_FILE"
  fi
  printf '[%s] Ejecutable fijado: SHA-256 %s.\n' "$(timestamp)" "$campaign_hash"
}

select_thread_budget() {
  local available candidate value
  available="$(nproc)"
  for value in "${SLURM_CPUS_PER_TASK:-}" "${SLURM_CPUS_ON_NODE:-}" "${PBS_NP:-}"; do
    candidate="${value%%(*}"
    if [[ "$candidate" =~ ^[0-9]+$ ]] && (( candidate < available )); then
      available="$candidate"
    fi
  done

  TOTAL_THREADS=0
  [[ "$PRIMARY_THREADS" =~ ^[1-9][0-9]*$ ]] || die 'PRIMARY_THREADS debe ser entero positivo'
  [[ "$FALLBACK_THREADS" =~ ^[0-9]+$ ]] || die 'FALLBACK_THREADS debe ser entero no negativo'

  if (( available >= PRIMARY_THREADS )); then
    TOTAL_THREADS="$PRIMARY_THREADS"
  elif (( FALLBACK_THREADS > 0 && available >= FALLBACK_THREADS )); then
    TOTAL_THREADS="$FALLBACK_THREADS"
  else
    if (( FALLBACK_THREADS > 0 )); then
      die "se requieren $PRIMARY_THREADS hilos, o $FALLBACK_THREADS como respaldo; disponibles: $available"
    fi
    die "se requieren exactamente $PRIMARY_THREADS hilos; disponibles: $available"
  fi

  [[ "$CONCURRENT_CASES" =~ ^[1-9][0-9]*$ ]] || die 'CONCURRENT_CASES debe ser entero positivo'
  (( TOTAL_THREADS % CONCURRENT_CASES == 0 )) || \
    die "$TOTAL_THREADS hilos no son divisibles entre $CONCURRENT_CASES casos simultáneos"
  THREADS_PER_CASE=$((TOTAL_THREADS / CONCURRENT_CASES))
  printf '%s\n' "$TOTAL_THREADS" >"$STATE_ROOT/total_threads.txt"
  printf '%s\n' "$THREADS_PER_CASE" >"$STATE_ROOT/threads_per_case.txt"
  printf '[%s] Presupuesto: %d hilos totales; %d casos simultáneos; %d hilos por caso.\n' \
    "$(timestamp)" "$TOTAL_THREADS" "$CONCURRENT_CASES" "$THREADS_PER_CASE"
}

declare -a CPU_SETS=()
OMP_PLACES_VALUE=cores
AFFINITY_UNIT='núcleos físicos'

configure_cpu_affinity() {
  local cpu core socket online key cpu_set slot index allowed_list available_count
  local -a physical_cpus=() logical_cpus=() selected_cpus=()
  local -a sockets=() sorted_sockets=() set_cpus=()
  local -A topology_cores=() seen_cores=() seen_sockets=()
  local -A physical_socket_sets=() logical_socket_sets=()

  case "$CPU_BINDING" in
    none)
      printf '[%s] Afinidad explícita desactivada.\n' "$(timestamp)"
      return 0
      ;;
    physical|smt) ;;
    *) die "CPU_BINDING debe ser none, physical o smt; recibido: $CPU_BINDING" ;;
  esac

  command -v lscpu >/dev/null 2>&1 || die 'la afinidad de CPU requiere lscpu'
  command -v taskset >/dev/null 2>&1 || die 'la afinidad de CPU requiere taskset'

  while IFS=, read -r cpu core socket online; do
    [[ "$cpu" =~ ^[0-9]+$ && "$core" =~ ^[0-9]+$ && "$socket" =~ ^[0-9]+$ ]] || continue
    key="$socket:$core"
    topology_cores[$key]=1
    # taskset descarta tanto CPU fuera del cpuset actual como CPU fuera de línea;
    # no dependemos del texto localizado de la columna Online de lscpu.
    taskset -c "$cpu" true >/dev/null 2>&1 || continue

    logical_cpus+=("$cpu")
    if [[ -z "${seen_sockets[$socket]+x}" ]]; then
      seen_sockets[$socket]=1
      sockets+=("$socket")
      logical_socket_sets[$socket]="$cpu"
    else
      logical_socket_sets[$socket]+=",$cpu"
    fi

    if [[ -z "${seen_cores[$key]+x}" ]]; then
      seen_cores[$key]=1
      physical_cpus+=("$cpu")
      if [[ -z "${physical_socket_sets[$socket]+x}" ]]; then
        physical_socket_sets[$socket]="$cpu"
      else
        physical_socket_sets[$socket]+=",$cpu"
      fi
    fi
  done < <(LC_ALL=C lscpu -p=CPU,CORE,SOCKET,ONLINE)

  if [[ "$CPU_BINDING" == physical ]]; then
    selected_cpus=("${physical_cpus[@]}")
    OMP_PLACES_VALUE=cores
    AFFINITY_UNIT='núcleos físicos'
  else
    selected_cpus=("${logical_cpus[@]}")
    OMP_PLACES_VALUE=threads
    AFFINITY_UNIT='hilos lógicos SMT'
  fi
  available_count=${#selected_cpus[@]}

  if (( available_count < TOTAL_THREADS )); then
    allowed_list="$(awk '/^Cpus_allowed_list:/ {print $2}' /proc/self/status 2>/dev/null || true)"
    die "se solicitaron $TOTAL_THREADS $AFFINITY_UNIT; lscpu detectó ${#topology_cores[@]} núcleos físicos y taskset permite $available_count unidades dentro de Cpus_allowed_list=${allowed_list:-desconocida}"
  fi

  # Cuando coincide un proceso por socket, cada caso conserva su caché L3 local.
  if (( ${#sockets[@]} == CONCURRENT_CASES )); then
    mapfile -t sorted_sockets < <(printf '%s\n' "${sockets[@]}" | sort -n)
    for socket in "${sorted_sockets[@]}"; do
      if [[ "$CPU_BINDING" == physical ]]; then
        cpu_set="${physical_socket_sets[$socket]}"
      else
        cpu_set="${logical_socket_sets[$socket]}"
      fi
      IFS=, read -r -a set_cpus <<<"$cpu_set"
      (( ${#set_cpus[@]} == THREADS_PER_CASE )) || \
        die "el socket $socket ofrece ${#set_cpus[@]} $AFFINITY_UNIT; se requieren $THREADS_PER_CASE"
      CPU_SETS+=("$cpu_set")
    done
  else
    # Respaldo genérico: conjuntos disjuntos de las unidades seleccionadas.
    index=0
    for ((slot=0; slot<CONCURRENT_CASES; slot++)); do
      cpu_set=""
      for ((core=0; core<THREADS_PER_CASE; core++)); do
        [[ -z "$cpu_set" ]] || cpu_set+=,
        cpu_set+="${selected_cpus[$index]}"
        index=$((index + 1))
      done
      CPU_SETS+=("$cpu_set")
    done
  fi

  for ((slot=0; slot<CONCURRENT_CASES; slot++)); do
    printf '[%s] Afinidad slot %d: CPU %s (%d %s).\n' \
      "$(timestamp)" "$slot" "${CPU_SETS[$slot]}" "$THREADS_PER_CASE" "$AFFINITY_UNIT"
  done
}

declare -a CASES=()
for reconstruction_name in $RECONSTRUCTIONS; do
  case "$reconstruction_name" in
    weno3|mp5|weno5) ;;
    *) die "reconstrucción no admitida para la campaña: $reconstruction_name" ;;
  esac
  CASES+=(
    "ef_a00_${reconstruction_name}"
    "ks_a00_${reconstruction_name}"
    "ks_a02_${reconstruction_name}"
    "ks_a04_${reconstruction_name}"
    "ks_a06_${reconstruction_name}"
    "ks_a09_${reconstruction_name}"
  )
done
readonly CASES
(( ${#CASES[@]} > 0 )) || die 'RECONSTRUCTIONS no puede estar vacío'

case_properties() {
  local label="$1"
  CASE_RECONSTRUCTION="${label##*_}"
  CASE_SCHEME="${CASE_RECONSTRUCTION}_hlle"
  if [[ "$label" == ef_* ]]; then
    CASE_METRIC=eddington_finkelstein
    CASE_GEOMETRY=spherical
    CASE_SPIN=0.0
  else
    CASE_METRIC=kerr_schild
    CASE_GEOMETRY=spheroidal
    case "$label" in
      ks_a00_*) CASE_SPIN=0.0 ;;
      ks_a02_*) CASE_SPIN=0.2 ;;
      ks_a04_*) CASE_SPIN=0.4 ;;
      ks_a06_*) CASE_SPIN=0.6 ;;
      ks_a09_*) CASE_SPIN=0.9 ;;
      *) die "espín no reconocido en $label" ;;
    esac
  fi
}

write_cold_parameter_file() {
  local label="$1" parameter_file output_folder
  case_properties "$label"
  parameter_file="$PARAMETER_ROOT/$label.par"
  output_folder="results/$label/session_000"
  cat >"$parameter_file" <<EOF
# Campaña PPI de congreso: resolución 800x1x400, t_final=5000 M.
problem = fishbone_equatorial
restart = false

metric = $CASE_METRIC
geometry = $CASE_GEOMETRY
bh_mass = 1.0
spin = $CASE_SPIN
adiabatic_index = 1.333333333333333333

fm_inner_radius = 6.25
fm_pressure_max_radius = 9.25
fm_polytropic_constant = 0.0015

nx = 800
ny = 1
nz = 400
r_min = 1.2
r_max = 40.0
theta_min = 0.0
theta_max = 1.0
phi_min = 0.0
phi_max = 2.0
logarithmic_r = true

reconstruction = $CASE_RECONSTRUCTION
riemann_solver = hlle
tvd_limiter = mc
shock_sensor = false

final_time = 5000.0
cfl = 0.4
save_interval = 50.0
checkpoint_interval = 250.0

output_prefix = $label
output_folder = $output_folder
vtk_mapping = physical

extract_gw = true
extract_mdot = true
ppi_diagnostics = true
diagnostic_stride = 1000
extraction_stride = 1000
mass_monitor_stride = 0

perturbation = pressure_noise
perturbation_seed = 3435
perturbation_time = 1000.0
perturbation_amplitude = 0.01
perturbation_mode = 1.0
EOF
}

write_restart_parameter_file() {
  local label="$1" checkpoint="$2" session="$3" parameter_file
  parameter_file="$PARAMETER_ROOT/${label}_${session}.par"
  cat >"$parameter_file" <<EOF
restart = true
restart_file = $checkpoint
final_time = 5000.0
cfl = 0.4
save_interval = 50.0
checkpoint_interval = 250.0
output_prefix = $label
output_folder = results/$label/$session
vtk_mapping = physical
shock_sensor = false
extract_gw = true
extract_mdot = true
ppi_diagnostics = true
diagnostic_stride = 1000
extraction_stride = 1000
mass_monitor_stride = 0
perturbation = pressure_noise
perturbation_seed = 3435
perturbation_time = 1000.0
perturbation_amplitude = 0.01
perturbation_mode = 1.0
EOF
  printf '%s\n' "$parameter_file"
}

validate_parameter_files() {
  local label failures=0 log_file
  for label in "${CASES[@]}"; do
    write_cold_parameter_file "$label"
    log_file="$LOG_ROOT/${label}.check.log"
    if (cd "$CAMPAIGN_ROOT" && "$CAMPAIGN_EXECUTABLE" --check "parameters/$label.par") \
      >"$log_file" 2>&1; then
      printf '[%s] PAR PASS %s\n' "$(timestamp)" "$label"
    else
      printf '[%s] PAR FAIL %s; consulte %s\n' "$(timestamp)" "$label" "$log_file"
      failures=$((failures + 1))
    fi
  done
  (( failures == 0 )) || die "$failures archivos .par no pasaron --check"
}

latest_checkpoint() {
  local label="$1"
  find "$RESULT_ROOT/$label" -type f -name 'checkpoint_*.rst' -printf '%T@ %p\n' \
    2>/dev/null | sort -n | tail -n 1 | awk '{print $2}'
}

next_session_name() {
  local label="$1" last_name last_number
  last_name="$(find "$RESULT_ROOT/$label" -mindepth 1 -maxdepth 1 -type d \
    -name 'session_[0-9][0-9][0-9]' -printf '%f\n' 2>/dev/null | sort | tail -n 1)"
  if [[ -z "$last_name" ]]; then
    printf 'session_000\n'
  else
    last_number="${last_name#session_}"
    printf 'session_%03d\n' "$((10#$last_number + 1))"
  fi
}

merge_series() {
  local label="$1" scheme="$2" file_name="$3"
  local combined_dir combined_file work_file source_file first_time
  local -a source_files=()
  combined_dir="$RESULT_ROOT/$label/combined"
  combined_file="$combined_dir/$file_name"
  work_file="$combined_dir/.${file_name}.tmp"
  mkdir -p "$combined_dir"

  mapfile -t source_files < <(find "$RESULT_ROOT/$label" -mindepth 3 -maxdepth 3 \
    -type f -path "*/$scheme/$file_name" -print | sort)
  (( ${#source_files[@]} > 0 )) || return 0

  awk '/^#/ {print; exit}' "${source_files[0]}" >"$combined_file"
  for source_file in "${source_files[@]}"; do
    first_time="$(awk '!/^#/ && NF {print $1; exit}' "$source_file")"
    [[ -n "$first_time" ]] || continue
    awk -v cutoff="$first_time" '/^#/ {print; next} $1 < cutoff {print}' \
      "$combined_file" >"$work_file"
    awk '!/^#/ && NF {print}' "$source_file" >>"$work_file"
    mv -f "$work_file" "$combined_file"
  done
}

merge_case_outputs() {
  local label="$1"
  case_properties "$label"
  merge_series "$label" "$CASE_SCHEME" ppi_modes.dat
  merge_series "$label" "$CASE_SCHEME" global_diagnostics.dat
  merge_series "$label" "$CASE_SCHEME" perturbation_events.dat
  merge_series "$label" "$CASE_SCHEME" GW_signal.dat
  merge_series "$label" "$CASE_SCHEME" m_dot.dat
}

audit_case() {
  local label="$1" diagnostics events
  diagnostics="$RESULT_ROOT/$label/combined/global_diagnostics.dat"
  events="$RESULT_ROOT/$label/combined/perturbation_events.dat"
  [[ -s "$diagnostics" && -s "$events" ]] || return 1
  awk '
    BEGIN {ok=1; rows=0; previous=-1.0e300}
    /^#/ {next}
    {
      rows++
      if ($1 < previous) ok=0
      previous=$1
      if ($5 <= 0.0 || $6 <= 0.0 || $9 < 0.0 || $9 >= 1.0 || $11 != 0) ok=0
      for (i=1; i<=11; i++) if ($i ~ /[Nn][Aa][Nn]|[Ii][Nn][Ff]/) ok=0
    }
    END {
      if (rows == 0 || previous < 4999.999999 || previous > 5000.000001) ok=0
      exit(ok ? 0 : 1)
    }
  ' "$diagnostics" || return 1
  awk '
    BEGIN {events=0; ok=1}
    /^#/ {next}
    NF {events++; if ($1 < 1000.0 || $2 != 1 || $3 != 3435) ok=0}
    END {exit(events == 1 && ok ? 0 : 1)}
  ' "$events"
}

run_case() {
  local label="$1" case_root checkpoint session parameter_file relative_parameter
  local slot="${2:-0}" cpu_set='' log_file exit_file status=0
  case_root="$RESULT_ROOT/$label"
  mkdir -p "$case_root"

  if [[ -s "$case_root/.complete" ]]; then
    printf '[%s] SKIP %s: ya está completo.\n' "$(timestamp)" "$label"
    return 0
  fi

  checkpoint="$(latest_checkpoint "$label")"
  if [[ -n "$checkpoint" ]]; then
    checkpoint="${checkpoint#"$CAMPAIGN_ROOT/"}"
    session="$(next_session_name "$label")"
    parameter_file="$(write_restart_parameter_file "$label" "$checkpoint" "$session")"
  else
    session=session_000
    parameter_file="$PARAMETER_ROOT/$label.par"
  fi
  relative_parameter="${parameter_file#"$CAMPAIGN_ROOT/"}"
  log_file="$LOG_ROOT/${label}_${session}.log"
  exit_file="$LOG_ROOT/${label}_${session}.exit"

  if ! (cd "$CAMPAIGN_ROOT" && "$CAMPAIGN_EXECUTABLE" --check "$relative_parameter") \
    >"$LOG_ROOT/${label}_${session}.check.log" 2>&1; then
    printf '2\n' >"$exit_file"
    printf '[%s] FAIL %s: el reinicio no pasó --check.\n' "$(timestamp)" "$label"
    return 1
  fi

  if [[ "$CPU_BINDING" != none ]]; then
    cpu_set="${CPU_SETS[$slot]}"
    printf '[%s] START %s, %s, %d %s, CPU %s, proceso inundaciones.\n' \
      "$(timestamp)" "$label" "$session" "$THREADS_PER_CASE" "$AFFINITY_UNIT" "$cpu_set"
  else
    printf '[%s] START %s, %s, %d hilos, proceso inundaciones.\n' \
      "$(timestamp)" "$label" "$session" "$THREADS_PER_CASE"
  fi
  printf '%s\n' "$(timestamp)" >"$LOG_ROOT/${label}_${session}.start"
  set +e
  if [[ "$CPU_BINDING" != none ]]; then
    (cd "$CAMPAIGN_ROOT" && \
      OMP_DYNAMIC=FALSE OMP_NUM_THREADS="$THREADS_PER_CASE" \
      OMP_PLACES="$OMP_PLACES_VALUE" OMP_PROC_BIND=CLOSE \
      taskset -c "$cpu_set" "$CAMPAIGN_EXECUTABLE" "$relative_parameter") \
      >"$log_file" 2>&1
  else
    (cd "$CAMPAIGN_ROOT" && \
      OMP_DYNAMIC=FALSE OMP_NUM_THREADS="$THREADS_PER_CASE" OMP_PROC_BIND=FALSE \
      "$CAMPAIGN_EXECUTABLE" "$relative_parameter") >"$log_file" 2>&1
  fi
  status=$?
  set -e
  printf '%s\n' "$status" >"$exit_file"
  printf '%s\n' "$(timestamp)" >"$LOG_ROOT/${label}_${session}.end"

  if (( status != 0 )); then
    printf '[%s] FAIL %s: salida %d; se reanudará al repetir el script.\n' \
      "$(timestamp)" "$label" "$status"
    return 1
  fi

  merge_case_outputs "$label"
  if ! audit_case "$label"; then
    printf '[%s] FAIL %s: terminó, pero no aprobó la auditoría final.\n' "$(timestamp)" "$label"
    return 1
  fi

  {
    printf 'case=%s\n' "$label"
    printf 'executable_sha256=%s\n' "$(awk '{print $1}' "$EXECUTABLE_HASH_FILE")"
    printf 'threads=%s\n' "$THREADS_PER_CASE"
    printf 'completed=%s\n' "$(timestamp)"
  } >"$case_root/.complete"
  printf '[%s] PASS %s\n' "$(timestamp)" "$label"
}

run_all_cases() {
  local label pid failures=0 slot=0
  local -a active_pids=() active_labels=()

  for label in "${CASES[@]}"; do
    run_case "$label" "$slot" &
    active_pids+=("$!")
    active_labels+=("$label")
    slot=$((slot + 1))

    if (( ${#active_pids[@]} == CONCURRENT_CASES )); then
      for pid in "${active_pids[@]}"; do
        if ! wait "$pid"; then failures=$((failures + 1)); fi
      done
      active_pids=()
      active_labels=()
      slot=0
    fi
  done

  for pid in "${active_pids[@]}"; do
    if ! wait "$pid"; then failures=$((failures + 1)); fi
  done
  return "$failures"
}

write_summary() {
  local label diagnostics last status
  printf 'case\tstatus\tt_final\tmass_total\tmass_disk\tangular_momentum\trho_min\tp_min\trho_max\tp_max\tv2_max\tn_atmosphere\tn_invalid\n' \
    >"$SUMMARY"
  for label in "${CASES[@]}"; do
    diagnostics="$RESULT_ROOT/$label/combined/global_diagnostics.dat"
    status=INCOMPLETE
    [[ -s "$RESULT_ROOT/$label/.complete" ]] && status=PASS
    if [[ -s "$diagnostics" ]]; then
      last="$(awk '!/^#/ && NF {line=$0} END {print line}' "$diagnostics")"
      printf '%s\t%s\t%s\n' "$label" "$status" "$(tr ' ' '\t' <<<"$last" | tr -s '\t')" \
        >>"$SUMMARY"
    else
      printf '%s\t%s\tNA\n' "$label" "$status" >>"$SUMMARY"
    fi
  done
}

printf '[%s] Inicio del controlador inundaciones.\n' "$(timestamp)"
prepare_executable
validate_parameter_files
select_thread_budget
configure_cpu_affinity

if [[ "$PREFLIGHT_ONLY" == 1 ]]; then
  printf '[%s] PREFLIGHT PASS: ejecutable, hilos y %d archivos .par.\n' \
    "$(timestamp)" "${#CASES[@]}"
  exit 0
fi

campaign_status=0
run_all_cases || campaign_status=$?
write_summary

if (( campaign_status != 0 )); then
  printf '[%s] Campaña incompleta: %d casos fallaron en esta sesión.\n' \
    "$(timestamp)" "$campaign_status"
  printf 'Repita el mismo comando para reanudar desde los checkpoints.\n'
  exit "$campaign_status"
fi

printf '[%s] CAMPAÑA COMPLETA: %d/%d PASS. Resumen: %s\n' \
  "$(timestamp)" "${#CASES[@]}" "${#CASES[@]}" "$SUMMARY"

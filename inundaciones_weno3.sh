#!/usr/bin/env bash
#SBATCH --job-name=inundaciones
#SBATCH --cpus-per-task=40
#SBATCH --hint=multithread
#PBS -N inundaciones
#PBS -l select=1:ncpus=40

# Perfil WENO3 de la campaña PPI. Reutiliza el controlador general probado.

set -Eeuo pipefail

readonly WRAPPER_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)/$(basename "${BASH_SOURCE[0]}")"
readonly WRAPPER_DIR="$(dirname "$WRAPPER_PATH")"
readonly WRAPPER_PARENT="$(cd "$WRAPPER_DIR/.." && pwd -P)"

export CAMPAIGN_ROOT="${CAMPAIGN_ROOT:-$WRAPPER_PARENT/inundaciones_ppi_weno3_800_t5000}"
export RECONSTRUCTIONS=weno3
export PRIMARY_THREADS="${PRIMARY_THREADS:-40}"
export FALLBACK_THREADS="${FALLBACK_THREADS:-0}"
export CONCURRENT_CASES="${CONCURRENT_CASES:-2}"
export CPU_BINDING="${CPU_BINDING:-smt}"
export LAUNCHER_NAME=inundaciones_weno3.sh
export CAMPAIGN_LABEL='PPI WENO3'

exec "$WRAPPER_DIR/inundaciones.sh" "$@"

#!/bin/bash -l
#SBATCH -A YOUR_PROJECT
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 01:00:00
#SBATCH -J MR_parallel
#SBATCH -o MR_parallel_%j.out
#SBATCH -e MR_parallel_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUR_EMAIL

set -u
set -o pipefail

module load bioinfo-tools augustus/3.3.3-CGP
module load python
module load R/4.3.1
module load R_packages/4.3.1
module load parallel/20210722-GCCcore-11.2.0

MR_SCRIPT="MR_NMR.R"
SUMMARY_SCRIPT="MR_results_arrange.py"

OUTPUT_DIR="./processed_files"
RESULTS_DIR="./results"
LOG_DIR="./logs"

mkdir -p "${OUTPUT_DIR}" "${RESULTS_DIR}" "${LOG_DIR}"
chmod 2775 "${OUTPUT_DIR}" "${RESULTS_DIR}" "${LOG_DIR}"

JOBS="${SLURM_NTASKS:-4}"

process_file() {
    local file_id="$1"
    local gwas_id="GCST${file_id}"
    local rel_path="${gwas_id}/harmonised/${gwas_id}.h.tsv.gz"
    local url="${BASE_URL}/${rel_path}"
    local downloaded="${OUTPUT_DIR}/${gwas_id}.h.tsv.gz"
    local tmp_file="${downloaded}.tmp"
    local out_dir="${RESULTS_DIR}/${gwas_id}"

    echo "[$(date '+%F %T')] ${gwas_id}: download started"

    rm -f "${tmp_file}"

    if ! wget -q --tries=3 --timeout=60 "${url}" -O "${tmp_file}"; then
        echo "[$(date '+%F %T')] ${gwas_id}: download failed"
        rm -f "${tmp_file}"
        return 1
    fi

    mv "${tmp_file}" "${downloaded}"

    echo "[$(date '+%F %T')] ${gwas_id}: MR started"

    mkdir -p "${out_dir}"

    if ! Rscript --vanilla "${MR_SCRIPT}" "${downloaded}" "${out_dir}"; then
        echo "[$(date '+%F %T')] ${gwas_id}: MR failed"
        rm -f "${downloaded}"
        return 1
    fi

    rm -f "${downloaded}"

    echo "[$(date '+%F %T')] ${gwas_id}: done"
    return 0
}

run_range() {
    local start_id="$1"
    local end_id="$2"
    local base_url="$3"
    local tag="$4"

    export BASE_URL="${base_url}"
    export OUTPUT_DIR RESULTS_DIR MR_SCRIPT
    export -f process_file

    echo "[$(date '+%F %T')] Range ${start_id}-${end_id} started"

    parallel \
        -j "${JOBS}" \
        --delay 0.5 \
        --joblog "${LOG_DIR}/parallel_${tag}.log" \
        process_file ::: $(seq "${start_id}" "${end_id}")

    local status=$?

    if [ "${status}" -ne 0 ]; then
        echo "[$(date '+%F %T')] Range ${start_id}-${end_id} finished with failures"
    else
        echo "[$(date '+%F %T')] Range ${start_id}-${end_id} finished successfully"
    fi

    return "${status}"
}

overall_status=0

run_range \
    90301941 \
    90302000 \
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90301001-GCST90302000" \
    "GCST90301941_90302000" || overall_status=1

run_range \
    90302001 \
    90302173 \
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90302001-GCST90303000" \
    "GCST90302001_90302173" || overall_status=1

echo "[$(date '+%F %T')] Running result summary"

if ! python "${SUMMARY_SCRIPT}"; then
    echo "[$(date '+%F %T')] Result summary failed"
    overall_status=1
fi

if [ "${overall_status}" -ne 0 ]; then
    echo "[$(date '+%F %T')] Finished with errors. Check ${LOG_DIR}/parallel_*.log and SLURM stderr."
    exit 1
fi

echo "[$(date '+%F %T')] All jobs completed"
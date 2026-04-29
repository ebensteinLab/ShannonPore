#!/usr/bin/env bash
# nanoentropy v4 — populate REFERENCE_DIR with FASTAs + GTFs.
#
# Downloads UCSC reference FASTAs (hg38 / mm10) plus matching
# ncbiRefSeq GTFs, then runs `samtools faidx` on each FASTA to produce
# the .fai index nanoentropy needs.
#
# By default downloads ALL genomes nanoentropy supports (currently hg38
# and mm10). Pass `--genome hg38` or `--genome mm10` to limit.
#
# Usage:
#   bash scripts/setup_references.sh                # both hg38 + mm10
#   bash scripts/setup_references.sh --genome hg38  # only hg38
#   NANOENTROPY_REF_DIR=/big/disk/refs bash scripts/setup_references.sh
#
# Sizes (post-decompression, approximate):
#   hg38.fa  ~3.2 GiB  + hg38.ncbiRefSeq.gtf.gz ~40 MiB
#   mm10.fa  ~2.8 GiB  + mm10.ncbiRefSeq.gtf.gz ~27 MiB
#
# This script is idempotent: each download is skipped if the target
# file already exists with non-zero size. It will NOT re-download or
# clobber existing files.

set -euo pipefail

V4_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF_DIR="${NANOENTROPY_REF_DIR:-$V4_DIR/reference_files}"

GENOMES_TO_FETCH=("hg38" "mm10")
case "${1:-}" in
    --genome)
        case "${2:-}" in
            hg38|mm10) GENOMES_TO_FETCH=("$2") ;;
            *) echo "error: --genome must be hg38 or mm10" >&2; exit 2 ;;
        esac
        ;;
    -h|--help|help)
        sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//'
        exit 0
        ;;
esac

log()  { printf '\033[1;34m[refs]\033[0m %s\n' "$*" >&2; }
warn() { printf '\033[1;33m[warn]\033[0m %s\n' "$*" >&2; }
err()  { printf '\033[1;31m[err]\033[0m %s\n'  "$*" >&2; exit 1; }

command -v curl >/dev/null 2>&1     || err "curl not found"
command -v gunzip >/dev/null 2>&1   || err "gunzip not found"
command -v samtools >/dev/null 2>&1 || err "samtools not found (install nanoentropy env first)"

mkdir -p "$REF_DIR"
log "REFERENCE_DIR = $REF_DIR"

download_if_missing() {
    local url="$1" dest="$2"
    if [[ -s "$dest" ]]; then
        log "skip (exists): $(basename "$dest")"
        return
    fi
    log "downloading: $url"
    curl -L --fail --retry 3 -o "$dest.partial" "$url"
    mv "$dest.partial" "$dest"
}

for g in "${GENOMES_TO_FETCH[@]}"; do
    case "$g" in
        hg38)
            FA_GZ_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
            GTF_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz"
            FA="$REF_DIR/hg38.fa"
            GTF="$REF_DIR/hg38.ncbiRefSeq.gtf.gz"
            ;;
        mm10)
            FA_GZ_URL="https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz"
            GTF_URL="https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.ncbiRefSeq.gtf.gz"
            FA="$REF_DIR/mm10.fa"
            GTF="$REF_DIR/mm10.ncbiRefSeq.gtf.gz"
            ;;
    esac

    if [[ ! -s "$FA" ]]; then
        log "$g: fetching FASTA (compressed)..."
        download_if_missing "$FA_GZ_URL" "$FA.gz"
        log "$g: gunzip FASTA (this can take a minute)..."
        gunzip -f "$FA.gz"
    else
        log "$g: FASTA already in place: $FA"
    fi

    if [[ ! -s "$FA.fai" ]]; then
        log "$g: building FASTA index (samtools faidx)..."
        samtools faidx "$FA"
    fi

    download_if_missing "$GTF_URL" "$GTF"
done

log "done. Verify with: nanoentropy doctor"

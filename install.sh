#!/usr/bin/env bash
# nanoentropy v4 — idempotent installer
# Produces a bit-for-bit reproducible environment from pinned lockfiles.
#
# Strategy:
#   1. Ensure micromamba (no admin rights needed; falls back to mamba/conda).
#   2. Create env from conda-lock.yml if present, else environment.yml.
#   3. Install pip-only deps from requirements.txt as a belt-and-braces step.
#   4. pip install -e .  →  the `nanoentropy` console script lands on PATH.
#   5. Run `nanoentropy doctor` and `nanoentropy selftest` to confirm.

set -euo pipefail

V4_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="${NANOENTROPY_ENV_NAME:-nanoentropy_v4}"
LOCKFILE="$V4_DIR/conda-lock.yml"
ENV_YML="$V4_DIR/environment.yml"
REQ_TXT="$V4_DIR/requirements.txt"

log()  { printf '\033[1;34m[install]\033[0m %s\n' "$*" >&2; }
warn() { printf '\033[1;33m[warn]\033[0m %s\n' "$*" >&2; }
err()  { printf '\033[1;31m[err]\033[0m %s\n' "$*" >&2; exit 1; }

# ── Resolve a conda-family CLI ────────────────────────────────────────────
choose_conda_cli() {
    if   command -v micromamba >/dev/null 2>&1; then echo micromamba
    elif command -v mamba      >/dev/null 2>&1; then echo mamba
    elif command -v conda      >/dev/null 2>&1; then echo conda
    else
        log "Installing micromamba into ~/.local/bin (no admin rights needed)..."
        mkdir -p "$HOME/.local/bin"
        local arch os
        arch="$(uname -m)"
        case "$arch" in
            x86_64) arch=64 ;;
            aarch64|arm64) arch=aarch64 ;;
            *) err "Unsupported arch: $arch" ;;
        esac
        os="$(uname -s | tr '[:upper:]' '[:lower:]')"
        curl -L "https://micro.mamba.pm/api/micromamba/${os}-${arch}/latest" \
            | tar -xj -C "$HOME/.local" bin/micromamba
        export PATH="$HOME/.local/bin:$PATH"
        echo micromamba
    fi
}

CONDA_BIN="$(choose_conda_cli)"
log "Using conda CLI: $CONDA_BIN"

# ── Create env ────────────────────────────────────────────────────────────
if [[ -f "$LOCKFILE" ]]; then
    log "Creating env '$ENV_NAME' from conda-lock.yml..."
    if command -v conda-lock >/dev/null 2>&1; then
        conda-lock install --name "$ENV_NAME" "$LOCKFILE"
    else
        warn "conda-lock CLI not on PATH; trying pip install of conda-lock"
        "$CONDA_BIN" run -n base pip install --quiet conda-lock || true
        conda-lock install --name "$ENV_NAME" "$LOCKFILE" \
            || "$CONDA_BIN" env create -n "$ENV_NAME" -f "$ENV_YML"
    fi
elif [[ -f "$ENV_YML" ]]; then
    warn "conda-lock.yml missing; falling back to environment.yml (resolver may pick newer transitive deps)."
    log "Creating env '$ENV_NAME' from environment.yml..."
    "$CONDA_BIN" env create -n "$ENV_NAME" -f "$ENV_YML" || \
        "$CONDA_BIN" env update -n "$ENV_NAME" -f "$ENV_YML"
else
    err "Neither conda-lock.yml nor environment.yml present in $V4_DIR"
fi

RUN() { "$CONDA_BIN" run -n "$ENV_NAME" "$@"; }

# ── Pip belt-and-braces ───────────────────────────────────────────────────
if [[ -f "$REQ_TXT" ]]; then
    log "Installing/verifying pip deps from requirements.txt..."
    RUN pip install --no-deps -r "$REQ_TXT" \
        || warn "pip reported issues (some pkgs already installed via conda)"
fi

# ── Executable bits on shipped scripts ────────────────────────────────────
log "Ensuring executable bits on scripts..."
chmod +x "$V4_DIR/install.sh" 2>/dev/null || true
chmod +x "$V4_DIR/bin/nanoentropy" 2>/dev/null || true
chmod +x "$V4_DIR/clawteam/orchestrate_review.sh" 2>/dev/null || true
chmod +x "$V4_DIR/scripts/setup_references.sh" 2>/dev/null || true

# ── REFERENCE_DIR ─────────────────────────────────────────────────────────
REF_DIR="${NANOENTROPY_REF_DIR:-$V4_DIR/reference_files}"
log "Ensuring REFERENCE_DIR exists: $REF_DIR"
mkdir -p "$REF_DIR"

# Reference FASTAs are MULTI-GIGABYTE (hg38 ~3.2 GiB, mm10 ~2.8 GiB).
# Auto-downloading by default would hang most installs, fill disks, and is
# a really bad first-run experience — so the install is OPT-IN.
#
#   bash install.sh                                             # no FASTA download (default)
#   NANOENTROPY_DOWNLOAD_REFERENCES=1 bash install.sh           # download both genomes
#   NANOENTROPY_DOWNLOAD_REFERENCES=hg38 bash install.sh        # download just one
#   NANOENTROPY_REF_DIR=/path/to/your/refs bash install.sh      # use existing refs
HAS_HG38="$([[ -s "$REF_DIR/hg38.fa" && -s "$REF_DIR/hg38.fa.fai" ]] && echo yes || echo no)"
HAS_MM10="$([[ -s "$REF_DIR/mm10.fa" && -s "$REF_DIR/mm10.fa.fai" ]] && echo yes || echo no)"

DL_REQ="${NANOENTROPY_DOWNLOAD_REFERENCES:-0}"
if [[ "$HAS_HG38" == "yes" && "$HAS_MM10" == "yes" ]]; then
    log "Reference FASTAs already present (hg38 + mm10)."
elif [[ "$DL_REQ" == "0" || -z "$DL_REQ" ]]; then
    log "Reference FASTAs not downloaded (default)."
    log "  → To use the GUI / CLI you need a populated REFERENCE_DIR."
    log "  → Either set NANOENTROPY_REF_DIR to your existing refs, OR run:"
    log "      bash $V4_DIR/scripts/setup_references.sh                # both"
    log "      bash $V4_DIR/scripts/setup_references.sh --genome hg38  # one"
else
    case "$DL_REQ" in
        1|both|all)         DL_ARGS=() ;;
        hg38)               DL_ARGS=(--genome hg38) ;;
        mm10)               DL_ARGS=(--genome mm10) ;;
        *) warn "NANOENTROPY_DOWNLOAD_REFERENCES=$DL_REQ — expected 1|hg38|mm10; skipping download."
           DL_ARGS=("--noop") ;;
    esac
    if [[ "${DL_ARGS[0]:-}" != "--noop" ]]; then
        log "Downloading UCSC reference FASTAs + GTFs (multi-GiB; this can take 30+ min)..."
        if RUN bash "$V4_DIR/scripts/setup_references.sh" "${DL_ARGS[@]}"; then
            log "Reference download complete."
        else
            warn "Reference download incomplete. Re-run later:"
            warn "    bash $V4_DIR/scripts/setup_references.sh ${DL_ARGS[*]}"
        fi
    fi
fi
export NANOENTROPY_REF_DIR="$REF_DIR"

# ── Editable install (puts `nanoentropy` console script on PATH) ──────────
log "Installing nanoentropy package (editable)..."
RUN pip install --no-deps -e "$V4_DIR" \
    || warn "editable install failed; CLI still usable as 'python -m src.cli'"

# ── Bash completion ───────────────────────────────────────────────────────
COMPLETION_DIR="$HOME/.local/share/bash-completion/completions"
if [[ -f "$V4_DIR/bin/nanoentropy.bash-completion" ]]; then
    mkdir -p "$COMPLETION_DIR"
    cp "$V4_DIR/bin/nanoentropy.bash-completion" "$COMPLETION_DIR/nanoentropy"
    log "Bash completion installed → $COMPLETION_DIR/nanoentropy"
fi

# ── Doctor ────────────────────────────────────────────────────────────────
log "Running 'nanoentropy doctor'..."
if RUN python -m src.cli doctor; then
    log "Doctor passed."
else
    warn "Doctor reported issues. Review the table above; you may need to set NANOENTROPY_REF_DIR."
fi

# ── Selftest ──────────────────────────────────────────────────────────────
if [[ "${NANOENTROPY_SKIP_SELFTEST:-0}" != "1" ]]; then
    log "Running 'nanoentropy selftest' (synthetic BAM end-to-end)..."
    if RUN python -m src.cli selftest; then
        log "Selftest passed. Install is fully functional."
    else
        warn "Selftest reported issues. See the output above."
    fi
else
    log "Skipping selftest (NANOENTROPY_SKIP_SELFTEST=1)."
fi

log "Done."
log "  GUI:  ${CONDA_BIN} activate $ENV_NAME && streamlit run $V4_DIR/app.py"
log "  CLI:  ${CONDA_BIN} activate $ENV_NAME && nanoentropy --help"

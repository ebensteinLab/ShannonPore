#!/usr/bin/env bash
# shannonpore v4 — idempotent installer
# Produces a bit-for-bit reproducible environment from pinned lockfiles.
#
# Strategy:
#   1. Ensure micromamba (no admin rights needed; falls back to mamba/conda).
#   2. Create env from conda-lock.yml if present, else environment.yml.
#   3. Install pip-only deps from requirements.txt as a belt-and-braces step.
#   4. pip install -e .  →  the `shannonpore` console script lands on PATH.
#   5. Run `shannonpore doctor` and `shannonpore selftest` to confirm.

set -euo pipefail

V4_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="${SHANNONPORE_ENV_NAME:-shannonpore}"
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
        # --fail so a non-2xx response (e.g. CDN error page) aborts instead
        # of being piped into tar. Default pipefail behaviour catches a
        # broken curl now.
        curl -L --fail --retry 3 \
            "https://micro.mamba.pm/api/micromamba/${os}-${arch}/latest" \
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
    if ! "$CONDA_BIN" env create -n "$ENV_NAME" -f "$ENV_YML"; then
        log "Env exists or create failed; trying env update instead..."
        "$CONDA_BIN" env update -n "$ENV_NAME" -f "$ENV_YML" \
            || err "Failed to create or update env '$ENV_NAME'. See errors above."
    fi
else
    err "Neither conda-lock.yml nor environment.yml present in $V4_DIR"
fi

# Sanity-check that the env actually exists before any RUN call.
if ! "$CONDA_BIN" env list 2>/dev/null | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    err "Env '$ENV_NAME' was not created. Inspect the errors above and re-run."
fi

# PYTHONNOUSERSITE=1 prevents `~/.local/lib/python3.10/site-packages` from
# shadowing the conda env's packages — a common-and-confusing source of
# "module not found" / "wrong version" errors on shared HPC accounts.
RUN() { PYTHONNOUSERSITE=1 "$CONDA_BIN" run -n "$ENV_NAME" "$@"; }

# ── Pip belt-and-braces ───────────────────────────────────────────────────
if [[ -f "$REQ_TXT" ]]; then
    log "Installing/verifying pip deps from requirements.txt..."
    RUN pip install --no-deps -r "$REQ_TXT" \
        || warn "pip reported issues (some pkgs already installed via conda)"
fi

# ── Executable bits on shipped scripts ────────────────────────────────────
log "Ensuring executable bits on scripts..."
chmod +x "$V4_DIR/install.sh" 2>/dev/null || true
chmod +x "$V4_DIR/bin/shannonpore" 2>/dev/null || true
chmod +x "$V4_DIR/bin/shannonpore-gui" 2>/dev/null || true
chmod +x "$V4_DIR/scripts/setup_references.sh" 2>/dev/null || true

# ── REFERENCE_DIR ─────────────────────────────────────────────────────────
REF_DIR="${SHANNONPORE_REF_DIR:-$V4_DIR/reference_files}"
log "Ensuring REFERENCE_DIR exists: $REF_DIR"
mkdir -p "$REF_DIR"

# Reference FASTAs are MULTI-GIGABYTE (hg38 ~3.2 GiB, mm10 ~2.8 GiB).
# As of this version they are downloaded BY DEFAULT so the tool is
# usable immediately after install. Opt out with:
#
#   SHANNONPORE_SKIP_REFERENCES=1 bash install.sh           # skip entirely
#   SHANNONPORE_DOWNLOAD_REFERENCES=hg38 bash install.sh    # only hg38
#   SHANNONPORE_DOWNLOAD_REFERENCES=mm10 bash install.sh    # only mm10
#   SHANNONPORE_REF_DIR=/path/to/your/refs bash install.sh  # use existing
HAS_HG38="$([[ -s "$REF_DIR/hg38.fa" && -s "$REF_DIR/hg38.fa.fai" ]] && echo yes || echo no)"
HAS_MM10="$([[ -s "$REF_DIR/mm10.fa" && -s "$REF_DIR/mm10.fa.fai" ]] && echo yes || echo no)"

if [[ "$HAS_HG38" == "yes" && "$HAS_MM10" == "yes" ]]; then
    log "Reference FASTAs already present (hg38 + mm10)."
elif [[ "${SHANNONPORE_SKIP_REFERENCES:-0}" == "1" ]]; then
    warn "SHANNONPORE_SKIP_REFERENCES=1 — skipping FASTA download."
    warn "Populate manually later:"
    warn "    bash $V4_DIR/scripts/setup_references.sh"
else
    # Honour SHANNONPORE_DOWNLOAD_REFERENCES if user pinned a single
    # genome; otherwise default to BOTH (hg38 + mm10).
    DL_REQ="${SHANNONPORE_DOWNLOAD_REFERENCES:-1}"
    case "$DL_REQ" in
        1|both|all)  DL_ARGS=() ;;
        hg38)        DL_ARGS=(--genome hg38) ;;
        mm10)        DL_ARGS=(--genome mm10) ;;
        *)
            warn "SHANNONPORE_DOWNLOAD_REFERENCES=$DL_REQ — expected 1|hg38|mm10."
            warn "Defaulting to both."
            DL_ARGS=()
            ;;
    esac
    log ""
    log "Downloading UCSC reference FASTAs + GTFs (multi-GiB)."
    log "This can take 20–60 min depending on bandwidth. Skip with"
    log "  SHANNONPORE_SKIP_REFERENCES=1 bash install.sh"
    log "or grab just one genome with"
    log "  SHANNONPORE_DOWNLOAD_REFERENCES=hg38|mm10 bash install.sh"
    log ""
    if RUN bash "$V4_DIR/scripts/setup_references.sh" "${DL_ARGS[@]}"; then
        log "Reference download complete."
    else
        warn "Reference download incomplete. Re-run later:"
        warn "    bash $V4_DIR/scripts/setup_references.sh ${DL_ARGS[*]}"
    fi
fi
export SHANNONPORE_REF_DIR="$REF_DIR"

# ── Editable install (puts `shannonpore` console script on PATH) ──────────
log "Installing shannonpore package (editable)..."
RUN pip install --no-deps --no-user -e "$V4_DIR" \
    || warn "editable install failed; CLI still usable as 'python -m src.cli'"

# ── Bash completion ───────────────────────────────────────────────────────
COMPLETION_DIR="$HOME/.local/share/bash-completion/completions"
if [[ -f "$V4_DIR/bin/shannonpore.bash-completion" ]]; then
    mkdir -p "$COMPLETION_DIR"
    cp "$V4_DIR/bin/shannonpore.bash-completion" "$COMPLETION_DIR/shannonpore"
    log "Bash completion installed → $COMPLETION_DIR/shannonpore"
fi

# ── Symlink wrappers into a bin dir on PATH ──────────────────────────────
#
# The wrappers use `mamba run -n shannonpore …` so they don't need
# `mamba activate`. Symlinking them onto PATH makes `shannonpore` and
# `shannonpore-gui` globally callable.
#
# Path discovery: scan $PATH for a user-writable bin directory that is
# ALREADY on PATH. This is robust against network-mounted home setups
# where $HOME and the PATH-dir parent don't share a literal prefix.
# Falls back to $HOME/.local/bin (and warns) if nothing on PATH is
# user-writable.
USER_BIN=""
PATH_HAS_USER_BIN="no"

IFS=':' read -ra _PATH_DIRS <<< "$PATH"

# `pick_bin <pattern>` returns the first PATH entry that:
#   - matches the glob in $1   (e.g. "*/.local/bin")
#   - exists on disk
#   - is user-writable
#   - is not a system / tool-specific dir
# Skips Rust .cargo/bin, conda condabin, nvm, miniconda envs, lab tool dirs,
# /usr/*, /etc/*, /bin, /sbin, /opt/*, /snap/*.
pick_bin() {
    local glob="$1" d
    for d in "${_PATH_DIRS[@]}"; do
        [[ -z "$d" ]] && continue
        [[ -d "$d" && -w "$d" ]] || continue
        case "$d" in
            /usr/*|/etc/*|/bin|/sbin|/opt/*|/snap/*) continue ;;
            */.cargo/bin|*/condabin|*/.nvm/*|*/miniconda*/envs/*) continue ;;
        esac
        # shellcheck disable=SC2053
        [[ "$d" == $glob ]] && { echo "$d"; return 0; }
    done
    return 1
}

# Try in order of preference: most generic/expected first.
USER_BIN="$(pick_bin '*/.local/bin'      || true)"
[[ -z "$USER_BIN" ]] && USER_BIN="$(pick_bin '*/local/bin'  || true)"
[[ -z "$USER_BIN" ]] && USER_BIN="$(pick_bin '*/bin'        || true)"

if [[ -n "$USER_BIN" ]]; then
    PATH_HAS_USER_BIN="yes"
else
    # Nothing on PATH was suitable — fall back to ~/.local/bin and warn.
    USER_BIN="$HOME/.local/bin"
    mkdir -p "$USER_BIN"
fi

log "Symlinking wrappers into: $USER_BIN  (on PATH: $PATH_HAS_USER_BIN)"
for w in shannonpore shannonpore-gui; do
    src="$V4_DIR/bin/$w"
    dst="$USER_BIN/$w"
    [[ -L "$dst" || -e "$dst" ]] && rm -f "$dst"
    ln -s "$src" "$dst"
    log "  $w → $dst"
done

if [[ "$PATH_HAS_USER_BIN" != "yes" ]]; then
    warn ""
    warn "$USER_BIN is NOT on your PATH yet."
    warn "Add this line to your shell rc and open a new terminal:"
    warn "    export PATH=\"$USER_BIN:\$PATH\""
    warn "Or call the wrappers by absolute path:"
    warn "    $V4_DIR/bin/shannonpore-gui"
    warn ""
fi

# ── Initialise the user's shell so `mamba activate` works ────────────────
# mamba 2.x refuses to mutate the parent shell unless `shell init` has
# been run first; without this the install ends with a confusing
# "Shell not initialized" error the first time the user tries to
# `mamba activate shannonpore`. We init their shell rc file once,
# idempotently. The bin/shannonpore* wrappers don't need this — they
# use `mamba run` — but interactive `streamlit run app.py` after
# `mamba activate` does.
USER_SHELL="$(basename "${SHELL:-/bin/bash}")"
RC_FILE=""
case "$USER_SHELL" in
    bash) RC_FILE="$HOME/.bashrc" ;;
    zsh)  RC_FILE="$HOME/.zshrc"  ;;
    fish) RC_FILE="$HOME/.config/fish/config.fish" ;;
esac

if [[ -n "$RC_FILE" ]]; then
    if grep -qE "(>>> mamba initialize|mamba shell hook|mamba shell init|conda initialize)" \
            "$RC_FILE" 2>/dev/null; then
        log "Shell ($USER_SHELL) already initialised for $CONDA_BIN."
    else
        log "Initialising $USER_SHELL so 'mamba activate $ENV_NAME' works in new shells..."
        if "$CONDA_BIN" shell init --shell "$USER_SHELL" 2>&1 | tail -3; then
            log "  Done. Open a new terminal or run:"
            log "    eval \"\$($CONDA_BIN shell hook --shell $USER_SHELL)\""
        else
            warn "Auto-init failed. Run manually once:"
            warn "    $CONDA_BIN shell init --shell $USER_SHELL"
        fi
    fi
fi

# ── Doctor ────────────────────────────────────────────────────────────────
log "Running 'shannonpore doctor'..."
if RUN python -m src.cli doctor; then
    log "Doctor passed."
else
    warn "Doctor reported issues. Review the table above; you may need to set SHANNONPORE_REF_DIR."
fi

# ── Selftest ──────────────────────────────────────────────────────────────
if [[ "${SHANNONPORE_SKIP_SELFTEST:-0}" != "1" ]]; then
    log "Running 'shannonpore selftest' (synthetic BAM end-to-end)..."
    if RUN python -m src.cli selftest; then
        log "Selftest passed. Install is fully functional."
    else
        warn "Selftest reported issues. See the output above."
    fi
else
    log "Skipping selftest (SHANNONPORE_SKIP_SELFTEST=1)."
fi

log ""
log "═════════════════════════════════════════════════════════════════"
log " Install complete. From any shell, right now:"
log ""
log "   shannonpore-gui          ← launch the Streamlit GUI"
log "   shannonpore --help       ← the CLI"
log ""
log " Both work without 'mamba activate' because they call into the"
log " conda env via 'mamba run' under the hood."
log ""
if [[ "$PATH_HAS_USER_BIN" != "yes" ]]; then
    log " (~/.local/bin isn't on your PATH — see the warning above and"
    log "  use the full path '$V4_DIR/bin/shannonpore-gui' meanwhile.)"
fi
log "═════════════════════════════════════════════════════════════════"

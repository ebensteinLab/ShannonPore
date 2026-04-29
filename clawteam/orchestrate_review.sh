#!/usr/bin/env bash
# ClawTeam adversarial review orchestrator for nanoentropy v4.
#
# Re-skinned for a Streamlit GUI debug + UX task (vs the upstream
# "research paper review" template). 5 agents:
#   1. code-qa       — walks every Streamlit callback; runs pytest+ruff
#   2. dataflow-qa   — traces displayed values back to source bedgraphs
#   3. ux-reviewer   — brutal Nielsen-heuristic critique
#   4. defender      — classifies each critique VALID/PARTIAL/INVALID
#   5. fixer         — implements fixes; re-runs tests; signals completion
#
# Two-layer loop:
#   • Inner: code-qa + dataflow-qa + ux-reviewer attack → defender rebuts →
#            fixer applies → re-review until ACCEPT (max MAX_ROUNDS).
#   • Outer: tear down team, spawn fresh team. Goal = first-round ACCEPT
#            from a fresh team (max MAX_SUBS).
#
# Usage:
#   bash orchestrate_review.sh BASE_TEAM_NAME [MAX_SUBS] [MAX_ROUNDS]
#   e.g. bash orchestrate_review.sh v4-ux-rev 3 5

set -euo pipefail

BASE_TEAM="${1:?usage: $0 BASE_TEAM_NAME [MAX_SUBS] [MAX_ROUNDS]}"
MAX_SUBS="${2:-3}"
MAX_ROUNDS="${3:-5}"

V4_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PROMPTS="$V4_DIR/clawteam/prompts"
LOG="$V4_DIR/clawteam/orchestrator_log.txt"
mkdir -p "$V4_DIR/clawteam"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOG"; }

require() {
    if ! command -v "$1" >/dev/null 2>&1; then
        log "ERROR: required command not found: $1"
        log "       install ClawTeam: pip install --user clawteam (or see"
        log "       https://github.com/HKUDS/ClawTeam)"
        exit 1
    fi
}
require clawteam
require tmux

spawn_team() {
    local TEAM="$1"
    local SUB="$2"

    log "Spawning team '$TEAM' for submission $SUB..."
    clawteam team spawn-team "$TEAM" -d "v4 GUI review submission $SUB" 2>/dev/null || true

    # 1. code-qa
    clawteam spawn -t "$TEAM" -n code-qa --agent-type code-qa --no-workspace \
        --task "$(cat "$PROMPTS/code-qa.md" \
                  | sed "s|@V4_DIR@|$V4_DIR|g; s|@TEAM@|$TEAM|g")"

    # 2. dataflow-qa (replaces hallucination-qa)
    clawteam spawn -t "$TEAM" -n dataflow-qa --agent-type hallucination-qa --no-workspace \
        --task "$(cat "$PROMPTS/dataflow-qa.md" \
                  | sed "s|@V4_DIR@|$V4_DIR|g; s|@TEAM@|$TEAM|g")"

    # 3. ux-reviewer (Reviewer 2)
    clawteam spawn -t "$TEAM" -n reviewer --agent-type adversarial-reviewer --no-workspace \
        --task "$(cat "$PROMPTS/ux-reviewer.md" \
                  | sed "s|@V4_DIR@|$V4_DIR|g; s|@TEAM@|$TEAM|g")"

    # 4. defender
    clawteam spawn -t "$TEAM" -n defender --agent-type research-defender --no-workspace \
        --task "$(cat "$PROMPTS/defender.md" \
                  | sed "s|@V4_DIR@|$V4_DIR|g; s|@TEAM@|$TEAM|g")"

    # 5. fixer
    clawteam spawn -t "$TEAM" -n fixer --agent-type code-fixer --no-workspace \
        --task "$(cat "$PROMPTS/fixer.md" \
                  | sed "s|@V4_DIR@|$V4_DIR|g; s|@TEAM@|$TEAM|g")"

    sleep 5
    for agent in code-qa dataflow-qa reviewer defender; do
        tmux send-keys -t "clawteam-$TEAM:$agent" \
            "Begin review now. Submission $SUB. ROUND 1. Write files prefixed ${TEAM}_ and suffixed _R1." Enter \
            2>/dev/null || true
    done
}

run_inner_loop() {
    local TEAM="$1"
    local SUB="$2"

    for round in $(seq 1 "$MAX_ROUNDS"); do
        log "  [Sub $SUB] Round $round / $MAX_ROUNDS"

        if [[ "$round" -gt 1 ]]; then
            log "  [Sub $SUB] Fixer round $round..."
            tmux send-keys -t "clawteam-$TEAM:fixer" \
                "ROUND $round: implement fixes; touch $V4_DIR/clawteam/${TEAM}_FIXES_R${round}_DONE when finished." Enter \
                2>/dev/null || true
            local waited=0
            while [[ ! -f "$V4_DIR/clawteam/${TEAM}_FIXES_R${round}_DONE" && "$waited" -lt 1800 ]]; do
                sleep 30; waited=$((waited + 30))
            done
            if [[ ! -f "$V4_DIR/clawteam/${TEAM}_FIXES_R${round}_DONE" ]]; then
                log "  [Sub $SUB] Fixer timeout; unblocking."
                touch "$V4_DIR/clawteam/${TEAM}_FIXES_R${round}_DONE"
            fi
        fi

        for agent in code-qa dataflow-qa reviewer defender; do
            tmux send-keys -t "clawteam-$TEAM:$agent" \
                "ROUND $round: re-review. Write ${TEAM}_${agent^^}_R${round}.md." Enter \
                2>/dev/null || true
        done

        local waited=0
        while [[ ! -f "$V4_DIR/clawteam/${TEAM}_REVIEWER_R${round}.md" && "$waited" -lt 1200 ]]; do
            sleep 30; waited=$((waited + 30))
        done
        [[ ! -f "$V4_DIR/clawteam/${TEAM}_REVIEWER_R${round}.md" ]] && continue

        local decision
        decision=$(grep -ioE "ACCEPT|MINOR|MAJOR|REJECT" \
                   "$V4_DIR/clawteam/${TEAM}_REVIEWER_R${round}.md" | head -1)
        log "  [Sub $SUB] Round $round decision: $decision"
        if [[ "$decision" == "ACCEPT" ]]; then
            echo "$round" > "$V4_DIR/clawteam/${TEAM}_ACCEPTED_ROUND"
            return 0
        fi

        local rebuttal_wait=0
        while [[ ! -f "$V4_DIR/clawteam/${TEAM}_REBUTTAL_R${round}.md" && "$rebuttal_wait" -lt 600 ]]; do
            sleep 15; rebuttal_wait=$((rebuttal_wait + 15))
        done
    done
    return 1
}

cleanup_team() {
    local TEAM="$1"
    yes y | clawteam team cleanup "$TEAM" 2>/dev/null || true
}

# ────────────────────── Main outer loop ──────────────────────────────────
log "═════════════════════════════════════════════════════════════"
log " ClawTeam adversarial review for nanoentropy v4"
log " Target: FIRST-ROUND ACCEPT from a fresh team"
log " Max submissions: $MAX_SUBS  ·  Max rounds each: $MAX_ROUNDS"
log " v4 dir: $V4_DIR"
log "═════════════════════════════════════════════════════════════"

for sub in $(seq 1 "$MAX_SUBS"); do
    TEAM="${BASE_TEAM}-sub${sub}"
    log ""
    log "╔ SUBMISSION $sub / $MAX_SUBS ════════════════════════════"

    spawn_team "$TEAM" "$sub"
    if run_inner_loop "$TEAM" "$sub"; then
        accepted=$(cat "$V4_DIR/clawteam/${TEAM}_ACCEPTED_ROUND" 2>/dev/null || echo "?")
        log "Submission $sub: ACCEPTED on round $accepted"
        if [[ "$accepted" == "1" ]]; then
            log "═══ FIRST-ROUND ACCEPT — v4 is bulletproof ═══"
            cat > "$V4_DIR/clawteam/FINAL_OUTCOME.md" <<EOF
FIRST-ROUND ACCEPT on submission $sub ($(date))
Team: $TEAM
All review files prefixed: ${TEAM}_
EOF
            exit 0
        fi
        log "Accepted but not on round 1 — fixes baked in. Spawning fresh team."
    else
        log "Submission $sub: not accepted after $MAX_ROUNDS rounds. Spawning fresh team."
    fi
done

log "Pipeline complete: $MAX_SUBS submissions attempted; not first-round accepted."
cat > "$V4_DIR/clawteam/FINAL_OUTCOME.md" <<EOF
NOT FIRST-ROUND ACCEPTED after $MAX_SUBS submissions ($(date)).
Review all ${BASE_TEAM}-sub*_*.md files for remaining issues.
EOF
exit 1

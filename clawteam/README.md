# ClawTeam adversarial review for v4

This directory contains the orchestrator + agent prompts that drive a
ClawTeam adversarial review of nanoentropy v4. Re-skinned from the
upstream "research paper review" template (see
`~/.claude/templates/adversarial_review_team.md`) for a Streamlit GUI
debug + UX task.

## What this does

Five agents, two loops:

| Agent | Mission | Output |
|---|---|---|
| `code-qa` | Walk every Streamlit callback + pipeline. Run pytest/ruff/py_compile. Flag stale state keys, race conditions, swallowed exceptions, off-by-ones. | `${TEAM}_CODE_QA_R{N}.md` |
| `dataflow-qa` | Trace every displayed value back to its source bedgraph. Catches axis swaps, label/data mismatches, BED 0/1-based bugs. | `${TEAM}_DATAFLOW_R{N}.md` |
| `ux-reviewer` | Brutal Nielsen-heuristic critic. Verdict: ACCEPT/MINOR/MAJOR/REJECT. | `${TEAM}_REVIEWER_R{N}.md` |
| `defender` | Lead-engineer rebuttal. Classifies each critique VALID / PARTIAL / INVALID and lists fixes. | `${TEAM}_REBUTTAL_R{N}.md` |
| `fixer` | Implements every VALID/PARTIAL fix. Re-runs tests + streamlit smoke. | `${TEAM}_FIX_TRACKER.md` |

* **Inner loop** (max `MAX_ROUNDS`, default 5): code-qa + dataflow-qa +
  reviewer attack → defender rebuts → fixer implements → re-review until
  ACCEPT.
* **Outer loop** (max `MAX_SUBS`, default 3): tear down team, spawn a
  *fresh* team. Goal = a fresh team accepts on **Round 1**.

`FINAL_OUTCOME.md` is written when the pipeline finishes, indicating
whether the goal was met or summarising remaining issues.

## Run

```bash
# 1. Make sure ClawTeam is installed (https://github.com/HKUDS/ClawTeam)
which clawteam || pip install --user clawteam

# 2. From the v4 directory:
cd /home/EbensteinLab/NanoEbenstein/Uri/entropy/New_Entropy_Tool_GUI/v4
nohup bash clawteam/orchestrate_review.sh v4-rev 3 5 > clawteam/orchestrator.out 2>&1 &

# 3. Monitor:
tail -f clawteam/orchestrator_log.txt
tmux attach -t clawteam-v4-rev-sub1     # watch live
clawteam board live                     # ClawTeam dashboard

# 4. When done:
cat clawteam/FINAL_OUTCOME.md
```

## Files

```
clawteam/
├── orchestrate_review.sh       # main orchestrator (this script)
├── prompts/
│   ├── code-qa.md              # Streamlit-aware code-quality reviewer
│   ├── dataflow-qa.md          # GUI value ↔ source-data traceability
│   ├── ux-reviewer.md          # brutal Nielsen-heuristic critic
│   ├── defender.md             # corresponding-engineer rebuttal
│   └── fixer.md                # implements all VALID/PARTIAL fixes
└── README.md                   # this file
```

## Why fresh-eyes Round-1 acceptance?

A team that has fought through 5 review rounds may rubber-stamp the
result by inertia. Tearing down the team and spawning a fresh one tests
whether the codebase is actually bulletproof — i.e., whether a reviewer
with zero history accepts immediately. If they don't, real issues
remain.

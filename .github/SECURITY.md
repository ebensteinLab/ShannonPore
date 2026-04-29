# Security policy

shannonpore is a research tool maintained by a small team. We take security
reports seriously and appreciate responsible disclosure.

## Reporting a vulnerability

Please do **not** open a public GitHub issue for security problems.

Instead, email the details privately to:

- **uribertocchi@gmail.com**

Include, where possible:

- A description of the issue and its impact
- Steps to reproduce or a proof-of-concept
- The affected version (`shannonpore --version`) and environment
- Any suggested mitigations

## Supported versions

Security fixes are issued against the latest minor version on the `main`
branch. Older versions receive fixes only when feasible.

| Version | Supported          |
|---------|--------------------|
| 4.x     | Yes (active)       |
| < 4.0   | No                 |

## Response SLA

This project is maintained on a best-effort basis as part of academic
research. We aim to:

- Acknowledge a report within **7 days**
- Provide an initial assessment within **30 days**
- Ship a fix or coordinated disclosure plan as soon as practical

If your report is time-sensitive, please say so explicitly in the email.

## Scope

In scope:

- Code in this repository (`src/`, `app.py`, `tests/`, helper scripts)
- The packaged CLI (`shannonpore`) and the Streamlit GUI

Out of scope:

- Vulnerabilities in upstream dependencies (please report those upstream)
- Issues that require a hostile local user with write access to the install

Thank you for helping keep shannonpore and its users safe.

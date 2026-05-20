# RECTIFY Agent Startup

This repository has active multi-agent work across M1, H2, and Sherlock.
Before touching pipeline code, read these files in order:

1. `AGENT_FIXES.md` — active bug coordination log. Update it when you find or
   resolve a pipeline bug.
2. `HANDOFF.md` — current session state, open work, and cluster-specific next
   steps.
3. `CLAUDE.md` — repo-specific operating rules, architecture pointers, HPC
   constraints, and surgical staging rules.
4. `dev/reviews/wip_disposition_2026-05-20.md` — current dirty-worktree
   triage. Use it before commenting that the tree is dirty.

Hard rules:

- Do not use `git add -A` or `git add .`; stage explicit paths only.
- Do not commit Kevin's existing WIP unless explicitly asked or the disposition
  ledger classifies it as safe.
- For FASTQ on HPC, use chunking/split workflows. Do not run unchunked
  `rectify run-all <fastq>` on cluster filesystems.
- Keep M1/GitHub/H2/Sherlock checkouts synchronized after committed work.


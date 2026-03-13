These files are frozen copies of the current Julia outputs from `/analyses/results`.

Purpose:
- Preserve the existing Julia-produced survey uncertainty outputs as golden-reference regression targets.
- Give future R or Macebase Analysis integrations a fixed baseline for parity checks.

Refresh policy:
- Do not overwrite these files casually.
- If they must change, regenerate from the Julia workflow intentionally and update `SHA256SUMS` in the same commit with an explanation.

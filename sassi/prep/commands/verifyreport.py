"""VERIFYREPORT: run verification problems and write the verification report (SASSI-EDU extension).

``VERIFYREPORT,[sel],[file],[fig]`` is the console/GUI form of ``python -m sassi.verify.report``
(:mod:`sassi.verify.report`): it runs the selected verification problems (requirements section 6)
and writes a Markdown report with, per problem, the computed and reference values, the error, the
tolerance and pass/fail (requirements section 6.1 item 4), the module coverage matrix (section 6.5)
and the lead decisions of sections 7.16/7.17 that explain the informative comparisons.

* ``sel``: ``P0`` (default, as VERIFY), ``P1``, ``P2``, ``ALL``, ``FAST`` (every problem not
  registered as slow) or problem ids joined with ``+`` (``VP-01+VP-30``; ``1`` means ``VP-01``).
* ``file``: output file, default ``VERIFICATION_REPORT.md`` in the model directory (MDL) or the
  working directory.  The project's own manual ``docs/verification/VERIFICATION_MANUAL.md`` is
  produced by the module command line, not by this command.
* ``fig``: 1 also writes the summary figure (PNG, ``figures/vp_margins.png`` next to the report).

The command is an action: it changes no model data, so WRITE writes nothing for it.  Each problem
runs in its own temporary directory, deleted when the problem has finished (D-W2-11).
"""
from __future__ import annotations

from pathlib import Path

from ..registry import command


@command("VERIFYREPORT", max_args=3, tier="P0", cls="action")
def cmd_verifyreport(c):
    """VERIFYREPORT,[sel],[file],[fig]: run verification problems (sel = P0 default, P1, P2, ALL, FAST or
    ids joined by '+') and write the Markdown verification report (default VERIFICATION_REPORT.md)."""
    from ...verify import load_all
    from ...verify import report as rep

    sel = c.str(1, "P0").strip().upper() or "P0"
    fname = c.str(2, "").strip() or "VERIFICATION_REPORT.md"
    figures = c.int(3, default=0) == 1
    reg = load_all()
    only = tiers = None
    skip_slow = False
    if sel in ("P0", "P1", "P2"):
        tiers = [sel]
    elif sel == "FAST":
        skip_slow = True
    elif sel != "ALL":
        try:
            only = rep.resolve_ids(sel.replace(";", "+").split("+"), reg)
        except KeyError as exc:
            c.fail(f"{exc.args[0]} (VERIFY,LIST lists the problems)")
    out = Path(c.output_path(fname))
    n_done = [0]

    def progress(k, n, run):
        n_done[0] = n
        c.info(f"{run.id}: {run.status}" + (f" ({run.elapsed:.1f} s)" if not run.skipped else ""))

    path, runs = rep.generate_manual(out, only=only, tiers=tiers, skip_slow=skip_slow, figures=figures,
                                     progress=progress, command=f"VERIFYREPORT,{c.str(1, '')},{c.str(2, '')}",
                                     title="SASSI-EDU Verification Report", registry=reg)
    bad = [r.id for r in runs if r.status in (rep.STATUS_FAIL, rep.STATUS_ERROR)]
    msg = (f"VERIFYREPORT {sel}: {len(runs)} problems, {sum(1 for r in runs if r.status == rep.STATUS_PASS)} "
           f"passed, {len(bad)} not passing; report written to {path}")
    if bad:
        c.warn(msg + f" ({', '.join(bad)})")
    else:
        c.confirm(msg)
        c.info(msg)

"""VP-49: section-cut integration T-S1 / T-S2 (sassi/verify/problems/vp_cuts.py)."""
from __future__ import annotations

from sassi.verify import run_problem


def test_vp49_section_cut_integration(tmp_path):
    res = run_problem("VP-49", tmp_path)
    failed = [c for c in res.checks if not c.passed]
    msg = "\n".join(f"{c.quantity}: computed {c.computed:.6g}, reference {c.reference:.6g}, "
                    f"error {c.error:.3g} > {c.tolerance:.3g} {c.note}" for c in failed[:20])
    assert res.passed, msg + "\n" + "\n".join(res.notes)
    # every published T-S1 / T-S2 quantity is checked
    names = {c.quantity for c in res.checks}
    for q in ("T-S1 Area", "T-S1 Ixx", "T-S1 Iyy", "T-S1 [Szz +10/-10] My", "T-S1 [Szz +10/-10, Sxz 3] Fx",
              "T-S1 [Syz +1/-1] Mz", "T-S2 Area", "T-S2 Iyy", "T-S2 [Sy'y' +100/-100] My", "T-S2 [Sy'y' +100 both] Fz"):
        assert q in names, q

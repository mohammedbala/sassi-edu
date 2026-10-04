"""Re-runs the numerical checks quoted in ../R1_sassi_theory.md (section 9).
Usage:  python run_all.py      (needs numpy + scipy)
Each script prints its own comparison table."""
import runpy, os
here = os.path.dirname(os.path.abspath(__file__))
os.chdir(here)
for s in ['tlm.py', 'test2.py', 'test3.py', 'test4.py', 'test4b.py', 'test4c.py', 'test5.py',
          'test5b_interp_rank.py', 'test9_damping.py', 'test10_waas_kausel.py', 'test7.py', 'point2.py', 'point3.py', 'test6.py', 'test8.py']:
    print('\n' + '=' * 20, s, '=' * 20)
    runpy.run_path(s, run_name='__main__')

import matplotlib; matplotlib.use('Agg')
import warnings, time; warnings.filterwarnings('ignore')
import centriole_kinematic as ck
from centriole_kinematic import Geometry, solve

print(f"{'BOND_HARD':>10} {'k_bond':>7} {'time':>7} {'worst gap':>10} {'outer':>8}")
print("-"*48)
for hard, kb in ((0.0, 12), (5.0, 12), (20.0, 12), (5.0, 60), (20.0, 60)):
    ck.BOND_HARD = hard
    t = time.time()
    try:
        s = solve(Geometry(), k_bond=kb, max_nfev=4000, rupture=False)
        gap = max(v['gap'] for v in s.bond_force.values())
        print(f"{hard:10.1f} {kb:7d} {time.time()-t:6.1f}s {gap:9.4f}nm {s.outer_diameter:8.2f}", flush=True)
    except Exception as e:
        print(f"{hard:10.1f} {kb:7d}  failed: {type(e).__name__}", flush=True)

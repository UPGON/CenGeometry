import matplotlib; matplotlib.use('Agg')
import warnings; warnings.filterwarnings('ignore')
from centriole_kinematic import *
cases = [
 ("min symmetry",        dict(N_cw=5,  N_mt=5)),
 ("max symmetry",        dict(N_cw=15, N_mt=15)),
 ("extreme mismatch",    dict(N_cw=5,  N_mt=15)),
 ("reverse mismatch",    dict(N_cw=15, N_mt=5)),
 ("singlet",             dict(MTn=1)),
 ("tiny spoke",          dict(spoke_rod=5.0)),
 ("huge spoke",          dict(spoke_rod=200.0)),
 ("tiny base",           dict(base_length=5.0)),
 ("huge base",           dict(base_length=150.0)),
 ("min protofilaments",  dict()),
 ("max protofilaments",  dict()),
]
for name, kw in cases:
    try:
        g = Geometry(**kw)
        if name == "min protofilaments": g = set_param(g,'n_pf_A',8)
        if name == "max protofilaments": g = set_param(g,'n_pf_A',20)
        s = solve(g)
        m = summarise(s)
        print(f"  {name:20} diam {m['diameter_nm']:7.1f}  strain {m['joint_rms_deg']:6.2f}  "
              f"clash {m['n_clashes']:3d}  conv {m['converged']}")
    except Exception as e:
        print(f"  {name:20} *** {type(e).__name__}: {e}")

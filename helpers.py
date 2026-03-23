import numpy as np
from params import *

def classify_pocket(g, h, pfl):
    if g == 0 or g == n_pf:
        return SITE_SEAM

    pf_r = g
    pf_l = g - 1

    right_lo = int(h     < pfl[pf_r])
    right_hi = int(h + 1 < pfl[pf_r])
    left_lo  = int(h     < pfl[pf_l])
    left_hi  = int(h + 1 < pfl[pf_l])

    n_present = right_lo + right_hi + left_lo + left_hi

    if n_present == 4:
        return SITE_LATTICE

    elif n_present == 3:
        return SITE_EDGE3

    elif n_present == 2:
        # Two tubulins stacked vertically on the same PF (longitudinal edge)
        if right_lo and right_hi and not left_lo:
            return SITE_EDGELONG2
        elif left_lo and left_hi and not right_lo:
            return SITE_EDGELONG2
        # Two tubulins side by side across adjacent PFs (lateral/horizontal edge)
        elif right_lo and left_lo:
            return SITE_EDGELAT2
        elif right_hi and left_hi:
            return SITE_EDGELAT2

    elif n_present == 1:
        return SITE_SINGLE

    else:  # n_present == 0
        raise ValueError(f"classify_pocket called on empty pocket g={g}, h={h}")
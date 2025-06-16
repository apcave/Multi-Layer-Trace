import math

import sys, os
sys.path.append(os.path.abspath('./build'))

from composite_wave_response import makeComposite, properateWave, getRp, getTp, getRs, getTs


def get_rp():
    # Pressure waves can have the same frequency but different angles.
    # Sum them together for an approximate result to compare with Levesque.
    pr = getRp()
    p_out = complex(0.0, 0.0)
    for p in pr:
        p_out += complex(p[0], p[1])
    return p_out

def get_tp():
    pr = getTp()    
    p_out = complex(0.0, 0.0)
    for p in pr:
        p_out += complex(p[0], p[1])
    return p_out

def get_rs():
    pr = getRs()
    p_out = complex(0.0, 0.0)
    for p in pr:
        p_out += complex(p[0], p[1])
    return p_out

def get_ts():
    pr = getTs()
    p_out = complex(0.0, 0.0)
    for p in pr:
        p_out += complex(p[0], p[1])
    return p_out
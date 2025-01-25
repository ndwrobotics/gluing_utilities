from sage.all_cmdline import *

from compute_traces import frob_trace_functions, is_HperpH_reducible
from search_db import search_ec_by_traces, filter_reducible_case
from symplectic_test import pass_symplectic_test
import pickle

        
def curve_to_safe_string(C):
    filename = str(C)[len("Hyperelliptic Curve over Rational Field defined by "):]
    filename = filename.replace("^", "")
    filename = filename.replace("*", "")
    return filename

def elliptic_curve_to_safe_string(E):
    filename = str(E)[len("Elliptic Curve defined by "):-len(" over Rational Field")]
    filename = filename.replace("^", "")
    filename = filename.replace("*", "")
    return filename

def compatible_curves(C,N,l):
    curves = set()
    for func in frob_trace_functions(C,N,l,bound=300):
        for E in filter_reducible_case(C, N, l, list(search_ec_by_traces(func, N, 100, l))):
            curves.add(E)
    curves = list(curves)
    symplectic_result, symplectic_run = pass_symplectic_test(C,N,curves,l)
    filtered_curves = []
    for E, result in zip(curves, symplectic_result):
        if result:
            filtered_curves.append(E)
    print(f"Found {len(filtered_curves)} potentially compatible curves")
    print(f"for {l=} and C = " + curve_to_safe_string(C))
    return filtered_curves, symplectic_run

            
def pickle_compatible_curves(C,N,l,directory="."):
    curves, symplectic_run = compatible_curves(C,N,l)
    print(curves)
    filename = directory + "/" + curve_to_safe_string(C)
    filename += " | " + str(l) + ".curves"
    with open(filename, "wb") as f:
        pickle.dump({
            "curve": C,
            "conductor": N,
            "l": l,
            "symplectic" : symplectic_run,
            "gluable_curves": curves,
            "reducible": is_HperpH_reducible(C,N,l),
        }, f)
    
def unpickle_compatible_curves(C, l,directory="."):
    filename = directory + "/" + curve_to_safe_string(C)
    filename += " | " + str(l) + ".curves"
    with open(filename, "rb") as f:
        return pickle.load(f)
from sage.all_cmdline import * 
#sys.path.append(os.environ['HOME'] + "/lmfdb")
#from lmfdb import dbSAGE_EXTCODE = SAGE_ENV['SAGE_EXTCODE'] 

script_dir = os.path.dirname(os.path.abspath(__file__))

## Setup magma console for gluing
def magma_setup():
    print(script_dir)
    magma.chdir(script_dir)
    magma.attach_spec("../../CHIMP/CHIMP.spec")
    #magma.attach("symplectic_test.m")
    magma.attach("gluing_utils.m")
    magma.attach("fast_gluing.m")

    magma.eval(f"""
        SetVerbose("QuarticIso", 1);
        SetVerbose("QuarticRec", 1);
        SetVerbose("Gluing", 1);
        SetVerbose("CurveRec", 1);
        SetVerbose("EndoFind", 1);
    """)


# Given an elliptic curve X1, a genus two curve X2, and a prime l,
# find all possible gluings of X1 and X2 along l-torsion
# where we operate with precision `prec` in the gluing step
# and precision `hprec` in the step of finding H.
def glue(X1, X2, l, prec, hprec=None):
    if hprec is None:
        print("Precision for finding H: default")
        print("Precision for finding invariants:", prec)
        return magma.function_call("GlueFast", [X1, X2, l, prec])
    else:
        print("Precision for finding H:", hprec)
        print("Precision for finding invariants:", prec)
        return magma.function_call("GlueFast", [X1, X2, l, prec, hprec])
    
# Given an elliptic curve X1, a genus two curve X2, and a prime l,
# find all possible (Diximier-Ohno or Shioda) invariants
# of gluings of X1 and X2 along l-torsion
# where we operate with precision `prec` in the gluing step
# and precision `hprec` in the step of finding H.
def glue_invs(X1, X2, l, prec, hprec=None):
    if hprec is None:
        print("Precision for finding H: default")
        print("Precision for finding invariants:", prec)
        return magma.function_call("GlueFastInv", [X1, X2, l, prec])
    else:
        print("Precision for finding H:", hprec)
        print("Precision for finding invariants:", prec)
        return magma.function_call("GlueFastInv", [X1, X2, l, prec, hprec])

# Given an elliptic curve X1, a genus two curve X2,
# verify that X3 is a gluing of X1 and X2.
def verify_gluing(X1, X2, X3, prec=500):
    return magma.function_call("VerifyGluing", [X1, X2, X3, prec])

def test():
    print("UHHH")
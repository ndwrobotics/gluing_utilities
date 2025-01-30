from sage.all_cmdline import * 
#sys.path.append(os.environ['HOME'] + "/lmfdb")
#from lmfdb import dbSAGE_EXTCODE = SAGE_ENV['SAGE_EXTCODE'] 

script_dir = os.path.dirname(os.path.abspath(__file__))


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


#def magma_reset():
#    magma.console()
#    magma_setup()

#def glue_slow(X1, X2, l, prec):
    #Q = magma.function_call("RationalsExtra", [prec])
    #print(Q)
    #return magma.function_call("AllGeometricGluingsCC", [X1, X2, Q, l])
    #return magma.function_call("GlueSlow", [X1, X2, l, prec])
    
def glue(X1, X2, l, prec, hprec=None):
    if hprec is None:
        return magma.function_call("GlueFast", [X1, X2, l, prec])
    else:
        return magma.function_call("GlueFast", [X1, X2, l, prec, hprec])
        
def glue_invs(X1, X2, l, prec, hprec=None):
    if hprec is None:
        return magma.function_call("GlueFastInv", [X1, X2, l, prec])
    else:
        return magma.function_call("GlueFastInv", [X1, X2, l, prec, hprec])

def verify_gluing(X1, X2, X3, prec=500):
    return magma.function_call("VerifyGluing", [X1, X2, X3, prec])
from sage.all_cmdline import *   # import sage library

SAGE_EXTCODE = SAGE_ENV['SAGE_EXTCODE'] 

script_dir = os.path.dirname(os.path.abspath(__file__))
magma.chdir(script_dir)
magma.attach_spec("../../CHIMP/CHIMP.spec")
magma.attach("symplectic_test.m")


#input: curve C, conductor N, a list of elliptic curves Es, prime ell
#return: 
# 1. symplectic test results in a string of True and False.
# 2. whether the test was run for all curves or not
def pass_symplectic_test(C, N, Es, ell):
    if(ell == 2):
        return [True for E in Es], True
    if(ell >= 23): #too large
        return [True for E in Es], False
    if(len(Es) == 0):
        return [], True
    p = magma.function_call('FindSymplecticTestPrime', [C, N, Es[0], ell])
    print("Symplectic test prime: ", p)
    if(p == -Integer(1) ):
        return [True for E in Es], False
    test_result = magma.function_call('SymplecticTestBatch', [C, N, Es, ell])
    result = []
    good = True
    for x in test_result:
        if (ell % Integer(4)  == Integer(1)) == (x == Integer(1)):
            result.append(True)
        elif x == Integer(0):
            result.append(True)
            good = False
        else:
            result.append(False)
    return result, good
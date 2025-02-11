
//helper function to find a torsion point
intrinsic FindGoodTorsionPoint(E::CrvEll, l::RngIntElt, p::RngIntElt, F::FldFin) -> PtEll
    {Given elliptic curve E, prime l, prime p (which is characteristic of field F),
    find a point P in E[l](F) such that Frob(P) is not an integer multiple of P.}
    f<x> := DivisionPolynomial(ChangeRing(E, GF(p)), l);
    for factor in Factorization(f) do
        g := factor[1];
        if ((l-1) div 2) mod Degree(g) eq 0 then
            continue;
        end if;
        _, px := HasRoot(g, F);
        u, v := HyperellipticPolynomials(E);
        R<y> := PolynomialRing(F);
        _, py := HasRoot(y^2 + Evaluate(v,px)*y - Evaluate(u,px), F);
        return elt<E | px, py>;
    end for;
    return 0;
end intrinsic;


intrinsic FindSymplecticTestPrime(C::CrvHyp, N::RngIntElt, E::CrvEll, l::RngIntElt) -> RngIntElt
    {Given genus 2 curve C with conductor N, elliptic curve E, and a prime l,
    find a prime p suitable for symplectic test (i.e., the one that Frob_p
    acts on E[l] by matrix of the form [[1,a],[0,1]] for a != 0.)
    Returns -1 if such prime p is not found up to 10000.}
    NE := Conductor(E);
    for p in PrimesUpTo(10000) do
        if (N mod p eq 0) or (NE mod p eq 0) or (p eq l) then
            continue;
        end if;
        t := TraceOfFrobeniusDirect(E, p);
        if (t*t - 4*p) mod l eq 0 then
            E1 := ChangeRing(E, GF(p));
            C1 := ChangeRing(C, GF(p));
            f<x> := DivisionPolynomial(E1, l); 
            if (Modexp(ChangeRing(x, GF(p)), p^(l-1), f) ne x) then
                R<x> := PolynomialRing(GF(l));
                t1 := GF(l)!t;
                if ChangeRing(FrobeniusPoly(C1),GF(l)) mod (x-t/2)^4 ne 0 then
                    return p;
                end if;
            end if;
        end if;
    end for;
    return -1;
end intrinsic;


intrinsic SymplecticTest(C::CrvHyp, E::CrvEll, l::RngIntElt, p::RngIntElt) -> RngIntElt
    {Perform Symplectic Test on genus 2 curve C, elliptic curve E, prime l,
    and prime p (obtained by FindSymplecticTestPrime).}
    F<a> := GF(p^(l*(l-1)));
    E1 := ChangeRing(E,F);
    C1 := ChangeRing(C,F);
    J1 := Jacobian(C1);
    
    P := FindGoodTorsionPoint(E1,l,p,F);
    sigmaP := elt<E1 | P[1]^p, P[2]^p, P[3]^p>;
    w1 := WeilPairing(P, sigmaP, l);
    print "Elliptic Curve done";
    
    B := GetTorsionBasis(J1,l,1);
    w2 := 1;
    f := ChangeRing(FrobeniusPoly(ChangeRing(C,GF(p))), GF(l));
    g := ChangeRing(FrobeniusPoly(ChangeRing(E,GF(p))), GF(l));
    h := f div g;
    coeffs := Coefficients(ChangeRing(h, Integers()));
    
    
    while true do
        a := Random([0..l-1]);
        b := Random([0..l-1]);
        c := Random([0..l-1]);
        d := Random([0..l-1]);
        P := a*B[1] + b*B[2] + c*B[3] + d*B[4];
        frobP := Frobenius(P, GF(p));
        //#construct point that is invariant in the Jordan block part
        Q := coeffs[1]*P + coeffs[2]*frobP + coeffs[3]*Frobenius(frobP,GF(p));
        frobQ := Frobenius(Q, GF(p));
        w2 := WeilPairing(Q, frobQ, l);
        if w2 ne 1 then
            break;
        end if;
    end while;
    
    for i in [1..(l-1)/2] do
        if w1^(i^2) eq w2 then
            return 1;
        end if;
    end for;
    
    return -1;
end intrinsic;

intrinsic SymplecticTestBatch(C::CrvHyp, N::RngIntElt, Es::[CrvEll], l::RngIntElt) -> SeqEnum
    {Perform Symplectic Test on genus 2 curve C, a list of elliptic curves Es, prime l,
    and prime p (obtained by FindSymplecticTestPrime).}
    result := [];
    for E in Es do
        p := FindSymplecticTestPrime(C,N,E,l);
        if(p ne -1) then
            Append(~result, SymplecticTest(C,E,l,p));
        else
            Append(~result, 0);
        end if;
    end for;
    return result;
end intrinsic;
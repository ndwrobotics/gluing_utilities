intrinsic GlueSlow(X1::Crv, X2::Crv, ell::RngIntElt, prec::RngIntElt) -> .
    {Wrapper for finding all gluings with a particular precision}
    R := RationalsExtra(prec);
    return AllGeometricGluingsCC(BaseChange(X1, R), BaseChange(X2, R), R, ell);
end intrinsic;

intrinsic GlueFast(X1::Crv, X2::Crv, ell::RngIntElt, prec::RngIntElt) -> .
    {Wrapper for efficiently finding all gluings with a particular algebraization precision}
    R := RationalsExtra(prec);
    return AllGeometricGluingsCCEfficient(BaseChange(X1, R), BaseChange(X2, R), R, ell);
end intrinsic;

intrinsic GlueFast(X1::Crv, X2::Crv, ell::RngIntElt, prec::RngIntElt, hprec::RngIntElt) -> .
    {Wrapper for efficiently finding all gluings with a particular algebraization precision and a particular subgroup search precision}
    R := RationalsExtra(prec);
    return AllGeometricGluingsCCEfficient(BaseChange(X1, R), BaseChange(X2, R), R, ell : hprec:=hprec);
end intrinsic;

intrinsic GlueFastInv(X1::Crv, X2::Crv, ell::RngIntElt, prec::RngIntElt, hprec::RngIntElt) -> .
    {Wrapper for efficiently finding all gluings with a particular algebraization precision and a particular subgroup search precision}
    R := RationalsExtra(prec);
    return AllGeometricGluingsInvCCEfficient(BaseChange(X1, R), BaseChange(X2, R), R, ell : hprec:=hprec);
end intrinsic;

intrinsic GlueFastInv(X1::Crv, X2::Crv, ell::RngIntElt, prec::RngIntElt) -> .
    {Wrapper for efficiently finding all gluings with a particular algebraization precision}
    R := RationalsExtra(prec);
    return AllGeometricGluingsInvCCEfficient(BaseChange(X1, R), BaseChange(X2, R), R, ell);
end intrinsic;

intrinsic VerifyGluing(X1::Crv, X2::Crv, X3::Crv, prec::RngIntElt) -> BoolElt
    {Verifies that a claimed gluing is correct.}
    F := RationalsExtra(prec);
    
    X1 := ChangeRing(X1, F);
    X2 := ChangeRing(X2, F);
    X3 := ChangeRing(X3, F);
    print X3;
    HeuristicEndomorphismRepresentation(X3);
    SetVerbose("EndoFind", 0);
    HeuristicEndomorphismDescription(X3);
    foo := HeuristicDecompositionFactors(X3);
    
    if #foo lt 2 then
        print "No factors discovered. Consider increasing the precision.";
        return false;
    end if;
    if Genus(foo[1]) eq 1 then
        ReconstructedX1 := foo[1];
        ReconstructedX2 := foo[2];
    elif Genus(foo[2]) eq 1 then
        ReconstructedX1 := foo[2];
        ReconstructedX2 := foo[1];
    else
        print "Factors have wrong genera???";
        return false;
    end if;
    
    print "Reconstructed X1:";
    print ReconstructedX1;
    print "Reconstructed X2:";
    print ReconstructedX2;
    
    b, twist := IsQuadraticTwist(ReconstructedX1, X1);
    if not b then
        print "Reconstructed X1 is not isogenous to original X1";
        return false;
    else
        print "Twist: ";
        print twist;
    end if;
    
    
    
    Hom23 := GeometricHomomorphismRepresentation(PeriodMatrix(X2), PeriodMatrix(X3), F);
    Hom13 := GeometricHomomorphismRepresentation(PeriodMatrix(X1), PeriodMatrix(X3), F);
    Hom32 := GeometricHomomorphismRepresentation(PeriodMatrix(X3), PeriodMatrix(X2), F);
    Hom31 := GeometricHomomorphismRepresentation(PeriodMatrix(X3), PeriodMatrix(X1), F);
    if (#Hom23 eq 0) or (#Hom13 eq 0) or (#Hom32 eq 0) or (#Hom31 eq 0) then
        print "Reconstructed factors are not actually factors.";
        return false;
    end if;
    print "Hom23", Hom23;
    print "Hom13", Hom13;
    print "Hom32", Hom32;
    print "Hom31", Hom31;
    return true;
end intrinsic;

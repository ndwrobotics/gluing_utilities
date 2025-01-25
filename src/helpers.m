intrinsic FindIdentityPrime(C::CrvHyp, l::RngIntElt, bound::RngIntElt) -> SeqEnum
    {Find a prime which Frob acts as identity matrix}
    D := Numerator(Discriminant(C));
    result := [];
    for p in PrimesUpTo(bound) do
        if (D mod p eq 0) or (p mod l ne 1) then
            continue;
        end if;
        B := GetTorsionBasis(Jacobian(ChangeRing(C,GF(p))),l,1);
        if (#B eq 4) then
            Append(~result, p);
        end if;
    end for;
    return result;
end intrinsic;

intrinsic IsFrobeniusMatrixDiagonal(C::CrvHyp, l::RngIntElt, p::RngIntElt) -> BoolElt
    {Genus 2}
    F<a> := GF(p^(l-1));
    C1 := ChangeRing(C,F);
    J1 := Jacobian(C1);
    B := GetTorsionBasis(J1,l,1);
    return (#B eq 4);
end intrinsic;

intrinsic IsFrobeniusMatrixDiagonal(E::CrvEll, l::RngIntElt, p::RngIntElt) -> BoolElt
    {Elliptic Curve}
    E1 := ChangeRing(E, GF(p));
    f<x> := DivisionPolynomial(E1, l);
    return Modexp(ChangeRing(x, GF(p)), p^(l-1), f) eq x;
end intrinsic;


intrinsic IsFrobeniusMatrixDiagonalBatch(C::CrvHyp, l::RngIntElt, primes::SeqEnum) -> SeqEnum
    {Genus 2 Curve. Return the result (1 if diagonal, 0 if not, -1 if invalid) for each prime p in primes.}
    result := [];
    d := Numerator(Discriminant(C));
    for p in primes do
        if d mod p eq 0 then
            Append(~result, -1);
        elif IsFrobeniusMatrixDiagonal(C, l, p) then
            Append(~result, 1);
        else
            Append(~result, 0);
        end if;
    end for;
    return result;
end intrinsic;

intrinsic IsFrobeniusMatrixDiagonalBatch(E::CrvEll, l::RngIntElt, primes::SeqEnum) -> BoolElt
    {Elliptic Curve. Return the result (1 if diagonal, 0 if not, -1 if invalid) for each prime p in primes.}
    result := [];
    d := Numerator(Discriminant(E));
    for p in primes do
        if d mod p eq 0 then
            Append(~result, -1);
        elif IsFrobeniusMatrixDiagonal(E, l, p) then
            Append(~result, 1);
        else
            Append(~result, 0);
        end if;
    end for;
    return result;
end intrinsic;
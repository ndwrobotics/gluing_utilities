/**
 * Code to efficiently construct all possible gluings of a genus 1 and genus 2 curve.
 * 
 * Modified from the gluing code written by Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 */


function LiftFF1(c, n)
    return (Integers() ! c)/n;
end function;

function QFromPVFor12(P, V)
/* Creates quotient of abelian variety corresponding to P by symplectic
 * subgroup corresponding to V */
p := Characteristic(BaseRing(V));
L1 := Lattice(IdentityMatrix(Rationals(), 6));
M2 := Matrix(Basis(V));
M2 := Matrix(Rationals(), #Rows(M2), #Rows(Transpose(M2)), [ LiftFF1(c, p) : c in Eltseq(M2) ]);
L2 := Lattice(M2);
L := L1 + L2;
T := Matrix(Basis(L));

E1 := StandardSymplecticMatrix(1);
E2 := StandardSymplecticMatrix(2);
E3 := DiagonalJoin(E1, E2);
T1 := Transpose(T);
X := Transpose(T1)*E3*T1;
E0, T2 := FrobeniusFormAlternating(ChangeRing(p*X, Integers()));
BT := T1*Transpose(T2);

Q := P*ChangeRing(BT, BaseRing(P));
assert IsBigPeriodMatrix(Q);
return Q;

end function;


// Identifies all Galois-stable one-dimensional subspaces of the ell-torsion of
// the Jacobian of the genus-2 curve C, with coordinates given w.r.t. the
// period matrix of the simplified model of Csimp
function OneDimSubspaces(C, ell : prec := 0)

oldprec := Precision(BaseField(C)`CC);
if prec eq 0 then
    prec := oldprec;
end if;
print "Precision for finding H:", prec;
CC := ComplexFieldExtra(prec);
R<T> := PolynomialRing(CC);
CSimp := SimplifiedModel(C);
A := AnalyticJacobian(CSimp : Precision := prec);
M := ChangeRing(PeriodMatrix(CSimp) * 2, CC);
P3 := ProjectiveSpace(GF(ell), 3);
generators := [];
for k -> Pt in Points(P3) do
    print "";
    print "---------------------------";
    print "H candidate number:", k, "/", #Points(P3);
    W := [];
    for i in [1..ell-1] do
        v := M * Matrix([ [CC! (Integers()!(i*Pt[j]))/ell ] : j in [1..4] ]);
        w := FromAnalyticJacobian(v, A);
        Append(~W, w);
    end for;
    p := R!1;
    for w in W do
        terms := [T - x[1] : x in w | Modulus(x[1]) lt 10^20 and Modulus(x[2]) lt 10^20];
        p := p * ((#terms eq 0) select 1 else (&*terms));
    end for;
    b,w := AlgebraizeElementsExtra(Coefficients(p), RationalsExtra(prec));
    if b then
        print "Found candidate H";
        Append(~generators, Pt);
    end if;
    print "---------------------------";
end for;
RationalsExtra(oldprec); // Necessary because Magma is bad and therefore RationalsExtra precision is global
return generators;

end function;


intrinsic InvariantsFromHCoordinates(P::ModMatFldElt, V::ModTupFld, VY::ModTupFld, a::Pt, F::Fld : index:=1, total:=1) -> .
{Takes as input a period matrix representing the product of a genus 1 and genus 2 curve, a 6-dimensional vector space V over F_ell, a 4-dimensional subspace VY of V, a point a in P3(F_ell), and a field F. Returns a list of invariants of potential ell-gluings of the two curves which use the coordinates of a as a generator for H and whose invariants are defined over the field F. (The optional parameters "index" and "total" are only used to control console output.)}
    FF := BaseField(V);
    h := V!([0,0] cat [a[i] : i in [1..4]]);

    H := sub<V | [h]>;

    Hperp := OrthogonalComplement(VY, H);
    HperpH, f := quo<Hperp|H>;

    h1 := (HperpH![1,0]) @@ f;
    h2 := (HperpH![0,1]) @@ f; //basis of Hperp
    h2 := h2 / (h1, h2); //make the basis symplectic
    Vs := [];
    for M in SL(2, FF) do
        u1 := h1 + V![M[1][1], M[1][2], 0, 0, 0, 0];
        u2 := h2 + V![-M[2][1], -M[2][2], 0, 0, 0, 0]; //#negative sign to make det = -1
        U := sub<V | [h,u1,u2]>;

        Append(~Vs, U);
    end for;
    invss := [];
    
    for i in [1..#Vs] do
        print "";
        print "-----------------------";
        print "Subgroup number:", i + (index - 1)*total, "/", total * #Vs;
        print "-----------------------";


        Q := QFromPVFor12(P, Vs[i]);
        tau := SmallPeriodMatrix(Q);

        invs, _, test := AlgebraizedInvariants(tau, F : Base := true);
        if test then
            print "";
            if #invs eq 9 then
                print "Shioda invariants found:";
            else
                print "Dixmier--Ohno invariants found:";
            end if;
            print invs;
            Append(~invss, invs);
        end if;
    end for;
    invss := [ invs : invs in Set(invss) ];
    return invss;
end intrinsic;

intrinsic AllGeometricGluingsCCEfficient(X1::Crv, X2::Crv, F::Fld, ell::RngIntElt : hprec:=300) -> .
{Finds all possible gluings of the genus 1 curve X1 and the genus 2 curve X2 along ell-torsion, with invariants defined over F. If the algorithm fails to identify H, hprec should be increased. If the algorithm failes to algebraize invariants, the precision of F should be increased.}

P1 := PeriodMatrix(SimplifiedModel(HyperellipticCurve(X1)));
P2 := PeriodMatrix(X2);
spaces := OneDimSubspaces(X2, ell : prec := hprec);

print "Coordinates for H candidates:";
print spaces;

P := DiagonalJoin(P1, P2);

F1 := StandardSymplecticMatrix(1);
F2 := StandardSymplecticMatrix(2);
F3 := DiagonalJoin(F1, F2);

FF := GF(ell);

V := VectorSpace(FF, 6, ChangeRing(F3, FF));
VX := sub<V | [V![1,0,0,0,0,0], V![0,1,0,0,0,0]]>;
VY := OrthogonalComplement(V, VX);

/* First find invariants */
invss := [ ];

Vs := [];

print "Precision for algebraizing invariants", Precision(F`CC);

for i -> a in spaces do
    invss cat:= InvariantsFromHCoordinates(P, V, VY, a, F : index:=i, total:=#spaces);
end for;

invss := [ invs : invs in Set(invss) ];

/* Then find curves */
Ys := [* *];
for invs in invss do
    if #invs eq 9 then
        Y := HyperellipticCurveFromShiodaInvariants(invs);
        //Y := ReducedMinimalWeierstrassModel(HyperellipticCurveFromShiodaInvariants(invs));
    else
        f, aut := TernaryQuarticFromDixmierOhnoInvariants(invs : minimize:=false);
        /*if #aut eq 0 then
            f := MinimizeC2Quartic(f);
        end if;*/
        Y := PlaneCurve(f);
    end if;
    Append(~Ys, Y);
end for;

return Ys;

end intrinsic;

intrinsic AllGeometricGluingsInvCCEfficient(X1::Crv, X2::Crv, F::Fld, ell::RngIntElt : hprec:=300) -> .
{Finds all possible invariatns of gluings of the genus 1 curve X1 and the genus 2 curve X2 along ell-torsion, with invariants defined over F. If the algorithm fails to identify H, hprec should be increased. If the algorithm failes to algebraize invariants, the precision of F should be increased.}

P1 := PeriodMatrix(SimplifiedModel(HyperellipticCurve(X1)));
P2 := PeriodMatrix(X2);
spaces := OneDimSubspaces(X2, ell : prec := hprec);

print "Coordinates for H candidates:";
print spaces;

P := DiagonalJoin(P1, P2);

F1 := StandardSymplecticMatrix(1);
F2 := StandardSymplecticMatrix(2);
F3 := DiagonalJoin(F1, F2);

FF := GF(ell);

V := VectorSpace(FF, 6, ChangeRing(F3, FF));
VX := sub<V | [V![1,0,0,0,0,0], V![0,1,0,0,0,0]]>;
VY := OrthogonalComplement(V, VX);

/* First find invariants */
invss := [ ];

Vs := [];

print "Precision for algebraizing invariants", Precision(F`CC);

for i -> a in spaces do
    invss cat:= InvariantsFromHCoordinates(P, V, VY, a, F : index:=i, total:=#spaces);
end for;

invss := [ invs : invs in Set(invss) ];

return invss;

end intrinsic;
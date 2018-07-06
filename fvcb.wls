#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Solve the FvCB model analytically
 *
 * Copyright (c) 2018 Wu Sun
 * License: MIT
 *)


(* An: net assimilation rate [\[Micro]mol m^-2 s^-1]
 * Ci: intercellular concentration of Subscript[CO, 2] [\[Micro]mol mol^-1]
 * Ca: ambient concentration of Subscript[CO, 2] [\[Micro]mol mol^-1]
 * gtot: total conductance of Subscript[CO, 2] [mol m^-2 s^-1]
 *)
An[Ci_] := gtot * (Ca - Ci)

(* Ac: carboxylation limited assimilation rate [\[Micro]mol m^-2 s^-1]
 * GammaStar: Subscript[CO, 2] compensation point [\[Micro]mol mol^-1]
 * Kco: effective Michaelis constant for RuBisCO [\[Micro]mol mol^-1]
 * Vcmax: maximum carboxylation rate [\[Micro]mol m^-2 s^-1]
 *)
Ac[Ci_] := (Ci - GammaStar) / (Ci + Kco) * Vcmax

(* Aj: electron transport limited assimilation rate [\[Micro]mol m^-2 s^-1]
 * J: actual electron transport rate [\[Micro]mol m^-2 s^-1]
 *)
Aj[Ci_] := (Ci - GammaStar) / (4 * Ci + 8 * GammaStar) * J

(* solve for the nontrivial intersection point of Ac and Aj *)
Cix = Ci /. Solve[Ac[Ci] == Aj[Ci] && Ci != GammaStar, Ci] // First

Print["(* Set global assumptions *)"]
globalAssum =
    Rd > 0 &&
    Ca > 0 &&
    Ci > 0 &&
    GammaStar > 0 &&
    Kco > 0 &&
    Vcmax > 0 &&
    J > 0 && J < 4 * Vcmax &&
    gtot > 0 &&
    Cix > 0
Print[globalAssum]

Print["\n(* Find conditions such that Cix > GammaStar *)"]
condCix = Reduce[{globalAssum, Cix > GammaStar}]
Print[condCix]

Print["\n(* Find conditions such that Aj > Ac *)"]
condCLimited = Reduce[{globalAssum, Cix > GammaStar, Aj[Ci] > Ac[Ci]}]
Print[condCLimited]

Print["\n(* Find conditions such that Aj < Ac *)"]
condELimited = Reduce[{globalAssum, Cix > GammaStar, Aj[Ci] < Ac[Ci]}]
Print[condELimited]

Print["\n(* CASE: carboxylation limited regime}\n(* Solve for Ci *)"]
Cic = Ci /. Solve[An[Ci] == Ac[Ci] - Rd, Ci]
Print["(* Root 1 *)\n", Cic[[1]], "\n(* Root 2 *)\n", Cic[[2]]]
Print["\n(* Calculate An results *)"]
Anc = Simplify[An[Cic]]
Print["(* Root 1 *)\n", Anc[[1]], "\n(* Root 2 *)\n", Anc[[2]]]

Print["\n(* CASE: electron transport limited regime *)\n(* Solve for Ci *)"]
Cij = Ci /. Solve[An[Ci] == Aj[Ci] - Rd, Ci]
Print["(* Root 1 *)\n", Cij[[1]], "\n(* Root 2 *)\n", Cij[[2]]]
Print["\n(* Calculate An results *)"]
Anj = Simplify[An[Cij]]
Print["(* Root 1 *)\n", Anj[[1]], "\n(* Root 2 *)\n", Anj[[2]]]

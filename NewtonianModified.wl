(* ::Package:: *)

BeginPackage["NewtonianVariableMass`"]


GravRates[ai_,bi_,ci_,af_,bf_,cf_,CM_]:=
Module[{ni0=ai,li0=bi,mi0=ci,nf0=af,lf0=bf,mf0=cf,M0=CM}, 
(rad=RadialInt[ni0,li0,nf0,lf0,M0]; sph=SphericalProd[li0,mi0,lf0,mf0]; 
Which[!Allowed[ni0,li0,mi0,nf0,lf0,mf0], "Invalid Rules",
mi0==mf0, SetPrecision[Power[TEnergy[ni0,nf0,M0],3]*Power[rad*GPZ[sph],2]/IntDem, sigfigs],
(mf0!=mi0), SetPrecision[Power[TEnergy[ni0,nf0,M0],3]*Power[Sqrt[Power[rad*GPX[sph],2]+Power[rad*GPY[sph],2]],2]/IntDem, sigfigs], 
True, "Default"])]; 


Begin["Private`"]


(* ::Text:: *)
(*Declare Private Functions*)
(**)
(*Reduced Mass*)


ReduceMass[a_,b_]:=a*b/(a+b);


Allowed[ni_,li_,mi_,nf_,lf_,mf_]:=Module[{ni0=ni,li0=li,mi0=mi,nf0=nf,lf0=lf,mf0=mf}, 
((ni0!=nf0) && (0<ni0) && (0<nf0) && (li0<ni0) && (lf0<nf0)  && (li0>= 0) && (lf0 >= 0)
 && (mi0 <= li0) && (mf0 <= lf0) && (mi0 >= -li0) && (mf0 >= -lf0) && 
(lf0-li0==1 || lf0-li0==-1) && (Abs[mf0-mi0]<= 1))];


(* ::Text:: *)
(*Energy of bound particle for <n>= a*)


ELevels[a_,M_]:=EFactor[M]/Power[a, 2];


(* ::Text:: *)
(*Transition Energies*)


TEnergy[a_,b_,M_]:=(ELevels[a,M]-ELevels[b,M])/\[HBar];


(* ::Text:: *)
(*Normalization constant <n>=a; <\[ScriptL]>=b*)


NormznConst[a_,b_,M_]:=Sqrt[Power[(2/(a*b0[M])), 3]*(a-b-1)!/(2a*Power[(a+b)!, 3])];


(* ::Text:: *)
(*Radial Component normalized <n>=a; <\[ScriptL]>=b*)


Rad[a_,b_,M_]:=NormznConst[a,b,M]*Exp[-r/(a*b0[M])]*Power[2r/(a*b0[M]), b]*(a+b)!*LaguerreL[a-b-1, 2b+1, 2r/(a*b0[M])]


(* ::Text:: *)
(*Radial Component Integrated  <ni>=ai; <\[ScriptL]i>=bi ; <nf>=af; <\[ScriptL]f>=bf*)


RadialInt[ai_,bi_,af_,bf_,M_]:=Integrate[echarge*Conjugate[Rad[ai,bi,M]]*Rad[af,bf,M]*Power[r,3], {r, 0, Infinity}];


(* ::Text:: *)
(* Spherical inner product  <\[ScriptL]i>=bi ; <mi>=ci; <\[ScriptL]f>=bf; <mf>=cf*)


SphericalProd[bi_,ci_,bf_,cf_]:=Conjugate[SphericalHarmonicY[bi,ci,\[Theta],\[Phi]]]*SphericalHarmonicY[bf,cf,\[Theta],\[Phi]]


(* ::Text:: *)
(*X, Y, Z Spherical Product Integrals*)


GPX[prodX_]:=Abs[Integrate[prodX*Power[Sin[\[Theta]], 2]*Cos[\[Phi]], {\[Phi], 0, 2Pi}, {\[Theta], 0, Pi}]]


GPY[prodY_]:=Abs[Integrate[prodY*Power[Sin[\[Theta]],2]*Sin[\[Phi]], {\[Phi], 0, 2Pi}, {\[Theta], 0, Pi}]]


GPZ[prodZ_]:=Abs[Integrate[prodZ*Sin[\[Theta]]*Cos[\[Theta]], {\[Phi], 0, 2Pi}, {\[Theta], 0, Pi}]]


(* ::Text:: *)
(*Reduced mass*)


rmass[M_]:=SetPrecision[ReduceMass[m, M],400]


(* ::Text:: *)
(*b0 parameter*)


b0[M_]:=SetPrecision[Power[\[HBar], 2]/(G*rmass[M]*m*M), 400];


(* ::Text:: *)
(*Energy level's constant factor*)


EFactor[M_]:=SetPrecision[-(rmass[M]*Power[G*m*M/\[HBar], 2])/2, 400];


(* ::Text:: *)
(*Declare Constants*)


(* ::Text:: *)
(*Number of significant figures*)


sigfigs=400;


(* ::Text:: *)
(*Electron mass in kilograms*)


(*ElectronMass=SetPrecision[9.10938188*10^-31, sigfigs];*) 


(* ::Text:: *)
(*Proton mass in kilograms*)


ProtonMass=SetPrecision[1.67262158*10^-27, sigfigs];


(* ::Text:: *)
(*Universal Gravitational constant Meters ^ 3 per (Kilograms * Seconds ^ 2) *)


G=SetPrecision[6.67398*10^-11, sigfigs];


(* ::Text:: *)
(*h bar Meters^2*Kilograms per Second*)


\[HBar]=SetPrecision[1.05457148*10^-34, sigfigs];


(* ::Text:: *)
(*Speed of Light (c) in Meters/Seconds*)


c=SetPrecision[299792458, sigfigs];


echarge=SetPrecision[1.60217657*10^-19, sigfigs]; 


(* ::Text:: *)
(*Permitivity of space in Farads/Meters*)


\[CurlyEpsilon]0=SetPrecision[8.854187817620*10^-12, sigfigs];


(* ::Text:: *)
(*Assign Variables*)


(* ::Text:: *)
(*Bound mass m in kilograms*)


m=ProtonMass;


(* ::Text:: *)
(*a0 Bohr radius*)


(*a0=SetPrecision[5.2917721092*10^-11, sigfigs];*)


(* ::Text:: *)
(*Overlap integral demonator constant factor*)


IntDem=SetPrecision[3*\[CurlyEpsilon]0*Pi*\[HBar]*Power[c, 3], sigfigs];


(*HEFactor=SetPrecision[-(rmass*Power[echarge, 4])/(32*Power[Pi*\[CurlyEpsilon]0*\[HBar], 2]), sigfigs];*)


End[]


EndPackage[]

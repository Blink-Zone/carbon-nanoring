(* ::Package:: *)

(* ::Chapter:: *)
(*NEGF*)


(* ::Section:: *)
(*Surface Green Function*)


(* ::Text:: *)
(*Eur. Phys. J. B 62, 381\[Dash]404 (2008).*)
(*I. Phys. F: Met. Phys. 15 (1985) 851-858.*)


(* ::Input::Initialization:: *)
\[Eta]=1.*^-7;
(*calculate isolated green function*)
Calcg0[\[Epsilon]_,h_]:=Inverse[(\[Epsilon]+I \[Eta])IdentityMatrix[Length[h]]-h];

(*iterate surface green function with Sancho method*)
IterSurG[{s_,e_,a_,b_},\[Epsilon]_]:=Module[{g},
g=Calcg0[\[Epsilon],e];
{
s+a . g . b,
e+a . g . b+b . g . a,
a . g . a,
b . g . b
}
]

(*calculate surface green function with Sancho method*)
CalcSurG[\[Epsilon]_,h00_,h01_,h10_]:=Module[{se,an},
se=NestWhile[IterSurG[#,\[Epsilon]]&,{h00,h00,h01,h10},(an=Norm[#[[3]]];
an>1.*^-8)&][[1]];
Calcg0[\[Epsilon],se]
]

(*iterate surface green function with dummy method*)
IterSurGD[g_,h00_,h01_,\[Epsilon]_]:=Calcg0[\[Epsilon],h00+h01\[ConjugateTranspose] . g . h01]
(*calculate surface green function with dummy method*)
CalcSurGD[\[Epsilon]_,h00_,h01_]:=Module[{g0,dg},
g0=Calcg0[\[Epsilon],h00];
NestWhile[IterSurGD[#,h01,h00,\[Epsilon]]&,g0,(dg=Abs[Im[Tr[#1-#2]]];(*PWL[dg];*)dg>1.*^-8)&,2,10000]
]

DOS[g_]:=-1/\[Pi] Im[Tr[g]]


(* ::Section:: *)
(*Surface Green Function with Overlapping*)


Calcgs0[\[Epsilon]_,h_,s_]:=Inverse[(\[Epsilon]+I \[Eta])s -h];
IterSurGS[{ss_,e_,a_,b_,s_},\[Epsilon]_]:=Module[{g},
g=Calcgs0[\[Epsilon],e,s];
{
ss+a . g . b,
e +a . g . b+b . g . a,
a . g . a,
b . g . b,
s
}
]

CalcSurGS[\[Epsilon]_,h00_,h01_,h10_,ol_]:=Module[{se,an},
se=NestWhile[IterSurGS[#,\[Epsilon]]&,{h00,h00,h01,h10,ol},(an=Norm[#[[3]]];an>1.*^-10)&][[1]];
Calcgs0[\[Epsilon],se,ol]
]


(*Dummy method*)
IterSurGSD[g_,h00_,h01_,ol_,\[Epsilon]_]:=Calcgs0[\[Epsilon],h00+h01\[ConjugateTranspose] . g . h01,ol]

CalcSurGSD[\[Epsilon]_,h00_,h01_,ol_]:=Module[{g0,dg},
g0=Calcg0[\[Epsilon],h00,ol];
NestWhile[IterSurGSD[#,h01,h00,ol,\[Epsilon]]&,g0,(dg=Abs[Im[Tr[#1-#2]]];(*PWL[dg];*)dg>1.*^-8)&,2,10000]
]


(* ::Section:: *)
(*Center Green Function*)


(* ::Input::Initialization:: *)
(*--------------------------NEGF--------------------------*)
(*EE[\[Epsilon]_]:=DiagonalMatrix[Table[\[Epsilon],{(NATOM-NELECL-NELECR)*4}]]*)
(*Old version abandoned*)
\[CapitalSigma]LK[\[Epsilon]_]:=VCL2 . gL[\[Epsilon]] . VLC2//SparseArray;
\[CapitalSigma]RK[\[Epsilon]_]:=VCR2 . gR[\[Epsilon]] . VRC2//SparseArray;
\[CapitalGamma]LK[\[Epsilon]_]:=I(\[CapitalSigma]LK[\[Epsilon]]-\[CapitalSigma]LK[\[Epsilon]]\[ConjugateTranspose]);
\[CapitalGamma]RK[\[Epsilon]_]:=I(\[CapitalSigma]RK[\[Epsilon]]-\[CapitalSigma]RK[\[Epsilon]]\[ConjugateTranspose]);
(*Green function Gr/a, dimension: 4N x 4N*)
(*Gr[\[Epsilon]_,\[Theta]_]:=Inverse[EE[\[Epsilon]]-HR[\[Theta]]-\[CapitalSigma]L-\[CapitalSigma]R];
Ga[\[Epsilon]_,\[Theta]_]:=Inverse[EE[\[Epsilon]]-HR[\[Theta]]+\[CapitalSigma]L+\[CapitalSigma]R];*)
Gr[\[Epsilon]_,h_]:=Calcg0[\[Epsilon],h+\[CapitalSigma]LK[\[Epsilon]]+\[CapitalSigma]RK[\[Epsilon]]];
Ga[\[Epsilon]_,h_]:=Gr[\[Epsilon],h]\[ConjugateTranspose];
(*spectra function AL/AR, dimension: 4N x 4N*)

Clear[AL,AR]
AL[\[Epsilon]_,h1_,h2_]:=(*AL[\[Epsilon],h1,h2]=*)Gr[\[Epsilon],h1] . \[CapitalGamma]LK[\[Epsilon]] . Ga[\[Epsilon],h2];
AR[\[Epsilon]_,h1_,h2_]:=(*AL[\[Epsilon],h1,h2]=*)Gr[\[Epsilon],h1] . \[CapitalGamma]RK[\[Epsilon]] . Ga[\[Epsilon],h2];

GC[\[Epsilon]_]:=Calcg0[\[Epsilon],hc+VCL2 . gL[\[Epsilon]] . VLC2+VCR2 . gR[\[Epsilon]] . VRC2]

Trans[\[Epsilon]_,h_]:=Gr[\[Epsilon],h] . \[CapitalGamma]RK[\[Epsilon]] . Ga[\[Epsilon],h] . \[CapitalGamma]LK[\[Epsilon]]//Tr//Re
(*\[CapitalGamma]LK[\[Epsilon]].*)(*Gr[\[Epsilon],h].(\[CapitalGamma]RK[\[Epsilon]]^2).Ga[\[Epsilon],h]//Tr//Re*)

(*--------------------------NEGF--------------------------*)


(* ::Section:: *)
(*NEGF trans*)


NEGFSystemWB\[CapitalSigma][\[Epsilon]_,h_,\[CapitalSigma]L_,\[CapitalSigma]R_]:=Module[{\[CapitalGamma]L,\[CapitalGamma]R,gr,ga,dos,trans},
{\[CapitalGamma]L,\[CapitalGamma]R}=-2Im[{\[CapitalSigma]L,\[CapitalSigma]R}];
gr=Calcg0[\[Epsilon],h+\[CapitalSigma]L+\[CapitalSigma]R];
ga=gr\[ConjugateTranspose];
dos=-1/\[Pi] Tr[Im[gr]];
trans=Tr[\[CapitalGamma]L . gr . \[CapitalGamma]R . ga];
{trans,dos}
]


NEGFSystemWB[\[Epsilon]_,h_,\[CapitalGamma]L_,\[CapitalGamma]R_]:=Module[{\[CapitalSigma]L,\[CapitalSigma]R,gr,ga,dos,trans},
{\[CapitalSigma]L,\[CapitalSigma]R}=-I/2 {\[CapitalGamma]L,\[CapitalGamma]R};

NEGFSystemWB\[CapitalSigma][\[Epsilon],h,\[CapitalSigma]L,\[CapitalSigma]R]
]


NEGFSystemET[\[Epsilon]_,HC_,HL_,HR_,tL_,tR_,tCL_,tCR_]:=Module[{GL,GR,\[CapitalSigma]L,\[CapitalSigma]R,\[CapitalGamma]L,\[CapitalGamma]R,gr,ga,dos,trans},
GL=CalcSurG[\[Epsilon],HL,tL,tL\[ConjugateTranspose]];
GR=CalcSurG[\[Epsilon],HR,tR,tR\[ConjugateTranspose]];
\[CapitalSigma]L=tCL . GL . tCL\[ConjugateTranspose];
\[CapitalSigma]R=tCR . GR . tCR\[ConjugateTranspose];

NEGFSystemWB\[CapitalSigma][\[Epsilon],HC,\[CapitalSigma]L,\[CapitalSigma]R]
]


(* ::Section:: *)
(*NEGF trans with overlap*)


NEGFSystemWB\[CapitalSigma]S[ee_,HC_,SC_,\[CapitalSigma]L_,\[CapitalSigma]R_]:=Module[{\[CapitalGamma]L,\[CapitalGamma]R,gr,ga,dos,trans},
\[CapitalGamma]L=-2Im[\[CapitalSigma]L];
\[CapitalGamma]R=-2Im[\[CapitalSigma]R];
gr=Calcgs0[ee,HC+\[CapitalSigma]L+\[CapitalSigma]R,SC];
ga=gr\[ConjugateTranspose];

dos=-1/\[Pi] Tr[Im[gr]];
trans=Tr[\[CapitalGamma]L . gr . \[CapitalGamma]R . ga];
{trans,dos}//Re
]


NEGFSystemETS[ee_,HC_,HL_,HR_,HCL_,HCR_,HL01_,HR01_,SC_,SL_,SR_,SCL_,SCR_,SL01_,SR01_]:=Module[
{tL,tR,tCL,tCR,GL,GR,\[CapitalSigma]L,\[CapitalSigma]R,\[CapitalGamma]L,\[CapitalGamma]R,gr,ga,dos,trans},
tL=HL01-ee SL01;
tR=HR01-ee SR01;
tCL=HCL-ee SCL;
tCR=HCR-ee SCR;

GL=CalcSurGS[ee,HL,tL,tL\[ConjugateTranspose],SL];
GR=CalcSurGS[ee,HR,tR,tR\[ConjugateTranspose],SR];
\[CapitalSigma]L=tCL . GL . tCL\[ConjugateTranspose];
\[CapitalSigma]R=tCR . GR . tCR\[ConjugateTranspose];

NEGFSystemWB\[CapitalSigma]S[ee,HC,SC,\[CapitalSigma]L,\[CapitalSigma]R]
]

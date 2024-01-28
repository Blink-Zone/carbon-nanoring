(* ::Package:: *)

(* ::Section:: *)
(*Message Control*)


On[Assert];


(* ::Section:: *)
(*Array Construct*)


BandArray[bands_,dims_]:=SparseArray[Band[#]->1&/@bands,dims];
ElementArray[elements_,dims_]:=SparseArray[#->1&/@elements,dims];
DiagSparseMatrix[list_,size_]:=DiagonalMatrix[Table[Boole[MemberQ[list,i]],{i,size}]];

ExpandI[m_,size_:4]:=Map[# IdentityMatrix[size]&,m,{2}]//ArrayFlatten;

BandArray2[row_]:=Module[{n},
(*right-down for first row*)
n=Length[row];
bands=Table[{Band[{1,i}]->row[[i]],Band[{i,1}]->row[[i]]},{i,n}]//Flatten;
SparseArray[bands,{n,n}]
]


(* ::Section:: *)
(*Derivate*)


(*\:77e9\:9635\:5217\:8868\:4e2d\:5fc3\:5dee\:5206\:6c42\:5bfc\:ff0c\:5faa\:73af\:5173\:8054\:4f1a\:6c42\:9996\:5c3e\:76f8\:8fde\:5012\:6570,\:8fb9\:754c\:533a\:57df\:6700\:597d\:4e3a\:96f6*)
MLCD[ml_,step_]:=Map[ListCorrelate[{-1/2,0,1/2},#,2]/step&,Transpose[ml],{2}]~Transpose~{2,3,1};


(* ::Section:: *)
(*Multiply*)


FastTrDot[a_,b_]:=Total[a*b\[Transpose],2]


(* ::Section:: *)
(*Hamiltonian Construct*)


GetH0DSDNA[nl_,{\[CapitalDelta]\[Epsilon]_,tin_,tic_}]:=Module[{},
h11=-\[CapitalDelta]\[Epsilon] IdentityMatrix[nl]+tin  BandArray[{{1,2},{2,1}},{nl,nl}];
h22=\[CapitalDelta]\[Epsilon] IdentityMatrix[nl]+tin BandArray[{{1,2},{2,1}},{nl,nl}];
h21=h12=tic IdentityMatrix[nl];
h={{h11,h12},{h21,h22}}//ArrayFlatten;
h
]

GetH0SSDNA[nl_,{\[Epsilon]0_,tin_}]:=Module[{},
h=\[Epsilon]0 IdentityMatrix[nl]+tin BandArray[{{1,2},{2,1}},{nl,nl}];
h
]


GetHelicalH2[atoms_,lc_,t1_]:=Module[{n},
(*Guo and Sun, PNAS 111,11658 (2014) more neighbor but without spin*)
n=Length[atoms];
(*neighbor distance dim:n-1 *)
li=Table[Norm[atoms[[i+1]]-atoms[[1]]],{i,3}];
ti=t1 Exp[-(li-li[[1]])/lc];
bands=Table[{Band[{1,1+i}]->ti[[i]],Band[{1+i,1}]->ti[[i]]},{i,n-1}]//Flatten;
SparseArray[bands,{n,n}]
]


(* ::Section:: *)
(*Electronic Structure*)


elementVE["H"]=1;
elementVE["C"]=4;
elementVE["N"]=5;
elementVE["O"]=6;
elementVE["S"]=6;
elementVE["Au"]=11;

GetHOLU[eval_,evec_,atoms_,types_]:=Module[{elementCount,NHOMO,NLUMO,EHOMO,ELUMO,EF,EGap},
(elementCount[#]=Count[types,#])&/@{"H","C","N","O","S","Au"};
NHOMO=Total[elementCount[#]*elementVE[#]&/@{"H","C","N","O","S","Au"}]/2//Ceiling;
NLUMO=NHOMO+1;
EHOMO=eval[[NHOMO]];
ELUMO=eval[[NLUMO]];
EF=(EHOMO+ELUMO)/2;
EGap=ELUMO-EHOMO;
{EHOMO,ELUMO,EF,EGap,NHOMO}
]


(* ::Section:: *)
(*Atom Coordinates Construct*)


GenHeliStruct[nrange_,la_,r_,\[Phi]_,\[Phi]0_]:=Module[{dz,dxy},
dxy=2 r Sin[\[Phi]/2];
dz=Sqrt[la^2-dxy^2];
Table[{r Cos[ i \[Phi]+\[Phi]0],r Sin[i \[Phi]+\[Phi]0], i dz},{i,nrange[[1]],nrange[[2]]}]//N
]


GenStraightStruct[nrange_,r0_,v_]:=Table[r0+i v,{i,nrange[[1]],nrange[[2]]}]//N


(* ::Section:: *)
(*Pauli*)


(*constants*)
\[Sigma]x={{0, 1},{1, 0}};
\[Sigma]y={{0, -I}, {I, 0}};
\[Sigma]z={{1, 0},{0, -1}};
\[Sigma]0={{0,0},{0,0}};


(* ::Section:: *)
(*PostProcess*)


PostProcessHATOM[h_,atoms_]:=Module[{},
Assert[Length[h]==Length[atoms]];
na=Length[atoms];
{x,y,z}=(atoms//Transpose);
{rx,ry,rz}=Outer[Subtract,#,#]&/@{x,y,z};
{ax,ay,az}=Outer[Plus,#,#]/2&/@{x,y,z};
{vx,vy,vz}=-I # h&/@{rx,ry,rz};
smoi=Sign[Im[ax vy- ay vx]];
lz=me(ax vy-ay vx);
]


(* ::Section:: *)
(*NEGF*)


Calcf[\[Epsilon]_,\[Mu]_,T_]:=Piecewise[{{1/(E^((\[Epsilon]-\[Mu])/T)+1),Abs[(\[Epsilon]-\[Mu])/T]<20}},Boole[\[Epsilon]<=\[Mu]]];
Calcdf[\[Epsilon]_,\[Mu]_,T_]:=Module[{f},f=Calcf[\[Epsilon],\[Mu],T];-1/T(f(1-f))]


Calcg0[ee_,h_]:=Inverse[(ee+I 1*^-8) IdentityMatrix[Length[h]]-h];


GetGE[h_,\[CapitalGamma]_,{emin_,emax_,estep_}]:=Module[{S,grE,gaE},
S=-I/2\[CapitalGamma];
grE=Table[Calcg0[ee,h+S],{ee,emin,emax,estep}];
gaE=#\[ConjugateTranspose]&/@grE;
{grE,gaE}
]


CorrelateGF[l1_,l2_]:=Module[{neh},
neh=(Length[l1]-1)/2;
ListCorrelate[l2,l1,neh+1,0][[neh+1;;-1]]
]
Get\[CapitalPi]lComp[{glE_,ggE_},{v1_,v2_},{i_,j_,k_,l_}]:=Module[{},
I v1[[l,k]]v2[[i,j]]CorrelateGF[glE[[All,k,i]],ggE[[All,j,l]]]//Re
]


(* ::Section:: *)
(*Photon Self*)


TrCorr[glList_,ggList_,m1_,m2_]:=Module[{n,neh,P1,P2,data,l1,l2,l3},
(*Tr[m1.gl.m2.gg]*)
n=Length[m1];
neh=(Length[glList]-1)/2;
P1=(m1 . #)&/@glList;
P2=(m2 . #)&/@ggList;
data=Table[
l1=P1[[All,i,j]];
l2=P2[[All,j,i]];
l3=ListCorrelate[l2,l1,neh+1,0][[neh+1;;-1]],
{i,n},{j,n}];
data~Total~2
];




(* ::Section:: *)
(*Plot Style*)


styles = {Frame -> True, AspectRatio -> 1, 
ImageSize -> {480, 480}, ImagePadding -> {{150, 30}, {120, 30}}, 
PlotStyle -> PointSize[Large], FillingStyle -> RGBColor[1, 0, 0], 
LabelStyle -> Directive[GrayLevel[0], 24, FontFamily -> "Times"], 
TicksStyle -> Directive[GrayLevel[0], 32, Thickness[Large], FontFamily -> "Times"], 
FrameStyle -> Directive[GrayLevel[0], 32, Thickness[Large], FontFamily -> "Times"]
};


(* ::Section:: *)
(*Eigen*)


MyEigensystem[m_]:=Module[{eigen},
eigen=Eigensystem[m]//N//Transpose;
eigen=SortBy[eigen,First]//Transpose;
(*Normalize*)
eigen[[2]]=Normalize/@ eigen[[2]];
eigen
]


(* ::Section:: *)
(*Plotting*)


PlotAtoms[atoms_] :=Module[{na},
    na = Length[atoms];
    Graphics3D[{Sphere /@ atoms, Table[Text[i, atoms[[i]]], {i, na}]}]
]


(* ::Section:: *)
(*Simultaneously Diagonalization*)


GetSDMatrix[a_,b_]:=Module[{p,q,bp,s},
p=a//Eigenvectors//Inverse;
bp=Inverse[p] . b . p//Chop;
q=bp//Eigenvectors//Inverse;
s=p . q;
s
]
(*GetSDMatrix2[a_,b_]:=Module[{p,q,bp,s},
p=Eigenvectors[a]\[ConjugateTranspose];
bp=p\[ConjugateTranspose].b.p;
q=Eigenvectors[bp]\[ConjugateTranspose];
s=p.q;
s
]*)


(* ::Section:: *)
(*Units*)


(* ::Text:: *)
(*\[HBar]=1,e=1,L=1\[Angstrom], E=1eV*)


(* ::Input::Initialization:: *)
unitT=Quantity[1,"ReducedPlanckConstant"/"Electronvolts"]//N//UnitConvert;
unitL=Quantity["Angstroms"]//N//UnitConvert;
unitM=Quantity[1,"Electronvolts"]/unitL^2*unitT^2//UnitConvert;
unitA=Quantity["ElementaryCharge"]/unitT//UnitConvert;
(*unitT=Quantity["Seconds"]//N//UnitConvert;
unitL=Quantity["Meters"]//N//UnitConvert;
unitM=Quantity["Kilograms"]//N//UnitConvert;
unitA=Quantity["Amperes"]//N//UnitConvert;*)
unitCollection={unitT,unitL,unitM,unitA}
GetMagnitude[q_]:=Module[{},
uList=UnitDimensions[q]//.{"ElectricCurrentUnit"->unitA,"LengthUnit"->unitL,"MassUnit"->unitM,"TimeUnit"->unitT};
unit=Times@@((#[[1]]^#[[2]])&/@uList);
q/unit
]
\[Epsilon]0=GetMagnitude[Quantity["ElectricConstant"]];
c=GetMagnitude[Quantity["SpeedOfLight"]];
\[HBar]=GetMagnitude[Quantity["ReducedPlanckConstant"]];
e=GetMagnitude[Quantity["ElementaryCharge"]];
Et=1;
q\[Alpha]=GetMagnitude[Quantity["FineStructureConstant"]]//N//UnitConvert;
{\[Epsilon]0,c,\[HBar],e,Et,q\[Alpha]}



With[{
qt=Quantity[1,"ReducedPlanckConstant"/"Electronvolts"],
 qeV=Quantity["Electronvolts"],
 qL=Quantity["Angstroms"], 
 qc=Quantity["SpeedOfLight"],
 q\[Epsilon]0=Quantity["ElectricConstant"],
 q\[HBar]=Quantity["ReducedPlanckConstant"],
 qe=Quantity["ElementaryCharge"],
 qme=Quantity["ElectronMass"]
 },
{
UnitConvert[q\[Epsilon]0,"Farads"/"Meters"],
GetMagnitude[qc],
GetMagnitude[q\[Epsilon]0],
UnitConvert[qeV^3  qL^2 qe^2/(3\[Pi]  q\[Epsilon]0 qc^3 q\[HBar]^4 )],
GetMagnitude[qe^2],
GetMagnitude[Quantity["Volts"]],
UnitConvert[qe^2 qeV^3 qL^2/(3\[Pi] q\[Epsilon]0 qc^3 q\[HBar]^4)],
GetMagnitude[qe^2 qeV^3 qL^2/(3\[Pi] q\[Epsilon]0 qc^3 q\[HBar]^4)],
GetMagnitude[qe^2/(3\[Pi] q\[Epsilon]0 qc^3)],
GetMagnitude[qe^2/(2\[Pi] q\[HBar])],
UnitConvert[qL/qt]//N,
UnitConvert[qe^2/(2\[Pi] q\[HBar]),"Siemens"]//N,
GetMagnitude[qme]
}
]



(* ::Section:: *)
(*Constant*)


\[Eta]=1*^-8;
me=Quantity["ElectronMass"]/(Quantity["Electronvolts"]/Quantity["Angstroms"]^2*(Quantity["ReducedPlanckConstant"/"Electronvolts"])^2);
me=GetMagnitude[Quantity["ElectronMass"]];

J0=2.1*^-8;


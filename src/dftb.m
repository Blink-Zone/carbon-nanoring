(* ::Package:: *)

(* ::Section:: *)
(*Parameters*)


(*\:5b9a\:4e49\:5404\:4e2a\:539f\:5b50\:7684\:6700\:5927\:8f68\:9053\:6570\:ff0c\:53ea\:80fd\:4e3a 1 4 9*)
MAXORB["H"]=1;
MAXORB["C"]=4;
MAXORB["N"]=4;
MAXORB["O"]=4;
MAXORB["S"]=4;
MAXORB["Au"]=9;
SPINFACTOR=1;
(*NEGF*)



(* ::Section:: *)
(*Units*)


(*Energy 1eV*)
uE=Quantity["Electronvolts"]//UnitConvert//N
(*Length 1\[Angstrom]*)
uL=Quantity["Angstroms"]//UnitConvert//N
(*Time h/eV*)
uT=Quantity["PlanckConstant"/"Electronvolts"]//UnitConvert//N
(*Mass Subscript[m, 0]*)
uM=Quantity["NeutronMass"]//UnitConvert//N


(* ::Section:: *)
(*XYZ File Manipulation*)


GetXYZTRAW[fileXYZ_]:=Module[{table,atoms,types},
table=Import[fileXYZ,"Table"][[3;;-1]];
types=table[[All,1]];
atoms=table[[All,2;;4]];
{atoms,types}
]
(*static Hamiltonian*)
(*HH:=Table[H[i,j],{i,NATOM},{j,NATOM}]~Flatten~{{1,3},{2,4}};*)

ImportXYZ[fileXYZ_]:=Module[{atoms,types,bonds},
{atoms,types}=GetXYZTRAW[fileXYZ];
bonds=GetBonds[atoms,types];
{atoms,types,bonds}
]

ImportMOL[fileMOL_]:=Module[{mol,atoms,types,bonds},
mol=Import[fileMOL];
atoms=mol["AtomCoordinates"]/Quantity["Angstroms"];
types=First/@AtomList[mol];
bonds=GetBonds[atoms,types];
{atoms,types,bonds}
]

ExportXYZ[path_,atoms_,types_]:=
Export[path,{types,atoms},{"XYZ",{"VertexTypes", "VertexCoordinates"}}]

SortXYZFile[f_]:=Module[{table,atoms,types},
table=Import[f,"Table"][[3;;-1]];
table=SortBy[table,{#[[4]]&,#[[2]]&,#[[3]]&}];
types=table[[All,1]];
atoms=table[[All,2;;4]];
CopyFile[f,f~~".bk"];
ExportXYZ[f,atoms,types]
]

MOL2XYZ[fileMOL_,fileXYZ_]:=Module[{atoms,types,bonds},
{atoms,types,bonds}=ImportMOL[fileMOL];
ExportXYZ[fileXYZ,atoms,types]
]


(* ::Section:: *)
(*DFTB*)


MAXBONDING["Au"]=4.0;
MAXBONDING["Others"]=4.0;

IsBond[atom1_,atom2_,type1_,type2_]:=\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
RowBox[{"0", "<", 
RowBox[{"Norm", "[", 
RowBox[{"atom1", "-", "atom2"}], "]"}], "<", 
RowBox[{"MAXBONDING", "[", "\"\<Au\>\"", "]"}]}], 
RowBox[{
RowBox[{"type1", "==", "\"\<Au\>\""}], "||", 
RowBox[{"type2", "==", "\"\<Au\>\""}]}]},
{
RowBox[{"0", "<", 
RowBox[{"Norm", "[", 
RowBox[{"atom1", "-", "atom2"}], "]"}], "<", 
RowBox[{"MAXBONDING", "[", "\"\<Others\>\"", "]"}]}], "True"}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\);

GetBonds[atoms_,types_]:=Module[{Natoms,array,bonds},
(*\:6ce8\:610f\:6b64\:51fd\:6570\:4ec5\:9002\:7528\:4e8eHCSAu\:4f53\:7cfb\:ff0c\:5047\:8bbe\:975e\:542b\:91d1\:952e\:952e\:957f<2A\:ff0c\:542b\:91d1\:952e\:952e\:957f<4A*)

Natoms=Length[atoms];
array=Table[{x,y},{y,Natoms},{x,Natoms}]~Flatten~1;
bonds=Select[array,IsBond[atoms[[#[[1]]]],atoms[[#[[2]]]],types[[#[[1]]]],types[[#[[2]]]]]&]
]




(* ::Section:: *)
(*Solve DFTB*)


(* ::Input::Initialization:: *)
(*DFTB\:6c42\:89e3\:5668*)
\[Xi]p=0.4;\[Xi]sp=0.6;

Slater[R_,a1_,a2_]:=Module[{bond=a1~~"-"~~a2,
bondR=a2~~"-"~~a1,
r=Normalize[R],
rr=Norm[R],
maxorb2=Max[MAXORB/@{a1,a2}],
l,m,n,ss\[Sigma],sp\[Sigma],sd\[Sigma],pp\[Sigma],pd\[Sigma],dd\[Sigma],pp\[Pi],pd\[Pi],dd\[Pi],dd\[Delta],ss\[Sigma]R,sp\[Sigma]R,sd\[Sigma]R,pp\[Sigma]R,pd\[Sigma]R,dd\[Sigma]R,pp\[Pi]R,pd\[Pi]R,dd\[Pi]R,dd\[Delta]R,pureHint},

{l,m,n}=r;

ss\[Sigma]=V[bond~~".ss\[Sigma]"][rr];
ss\[Sigma]R=V[bondR~~".ss\[Sigma]"][rr];
If[maxorb2>1,
sp\[Sigma]=V[bond~~".sp\[Sigma]"][rr];
pp\[Sigma]=V[bond~~".pp\[Sigma]"][rr];
pp\[Pi]=V[bond~~".pp\[Pi]"][rr];

sp\[Sigma]R=V[bondR~~".sp\[Sigma]"][rr];
pp\[Sigma]R=V[bondR~~".pp\[Sigma]"][rr];
pp\[Pi]R=V[bondR~~".pp\[Pi]"][rr];
];
If[maxorb2>4,
sd\[Sigma]=V[bond~~".sd\[Sigma]"][rr];
dd\[Sigma]=V[bond~~".dd\[Sigma]"][rr];
pd\[Sigma]=V[bond~~".pd\[Sigma]"][rr];
pd\[Pi]=V[bond~~".pd\[Pi]"][rr];
dd\[Pi]=V[bond~~".dd\[Pi]"][rr];
dd\[Delta]=V[bond~~".dd\[Delta]"][rr];
sd\[Sigma]R=V[bondR~~".sd\[Sigma]"][rr];
pd\[Sigma]R=V[bondR~~".pd\[Sigma]"][rr];
dd\[Sigma]R=V[bondR~~".dd\[Sigma]"][rr];
pd\[Pi]R=V[bondR~~".pd\[Pi]"][rr];
dd\[Pi]R=V[bondR~~".dd\[Pi]"][rr];
dd\[Delta]R=V[bondR~~".dd\[Delta]"][rr];
];

(*s x y z xy yz zx x^2-y^2 3z^2-r^2*)
If[maxorb2==1,pureHint={{ss\[Sigma]}},
If[maxorb2==4,pureHint=\!\(\*
TagBox[
TagBox[
TagBox[
RowBox[{"(", GridBox[{
{"ss\[Sigma]", 
RowBox[{"l", " ", "sp\[Sigma]"}], 
RowBox[{"m", " ", "sp\[Sigma]"}], 
RowBox[{"n", " ", "sp\[Sigma]"}]},
{
RowBox[{
RowBox[{"-", "l"}], " ", "sp\[Sigma]R"}], 
RowBox[{
RowBox[{
RowBox[{"(", 
RowBox[{"1", "-", 
SuperscriptBox["l", "2"]}], ")"}], " ", "pp\[Pi]"}], "+", 
RowBox[{
SuperscriptBox["l", "2"], " ", "pp\[Sigma]"}]}], 
RowBox[{"l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "pp\[Pi]"}], "+", "pp\[Sigma]"}], ")"}]}], 
RowBox[{"l", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "pp\[Pi]"}], "+", "pp\[Sigma]"}], ")"}]}]},
{
RowBox[{
RowBox[{"-", "m"}], " ", "sp\[Sigma]R"}], 
RowBox[{"l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "pp\[Pi]R"}], "+", "pp\[Sigma]R"}], ")"}]}], 
RowBox[{
RowBox[{
RowBox[{"(", 
RowBox[{"1", "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "pp\[Pi]"}], "+", 
RowBox[{
SuperscriptBox["m", "2"], " ", "pp\[Sigma]"}]}], 
RowBox[{"m", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "pp\[Pi]"}], "+", "pp\[Sigma]"}], ")"}]}]},
{
RowBox[{
RowBox[{"-", "n"}], " ", "sp\[Sigma]R"}], 
RowBox[{"l", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "pp\[Pi]R"}], "+", "pp\[Sigma]R"}], ")"}]}], 
RowBox[{"m", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "pp\[Pi]R"}], "+", "pp\[Sigma]R"}], ")"}]}], 
RowBox[{
RowBox[{
RowBox[{"(", 
RowBox[{"1", "-", 
SuperscriptBox["n", "2"]}], ")"}], " ", "pp\[Pi]"}], "+", 
RowBox[{
SuperscriptBox["n", "2"], " ", "pp\[Sigma]"}]}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\),
If[maxorb2==9,pureHint=\!\(\*
TagBox[
TagBox[
TagBox[
RowBox[{"(", GridBox[{
{"ss\[Sigma]", 
RowBox[{"l", " ", "sp\[Sigma]"}], 
RowBox[{"m", " ", "sp\[Sigma]"}], 
RowBox[{"n", " ", "sp\[Sigma]"}], 
RowBox[{
SqrtBox["3"], " ", "l", " ", "m", " ", "sd\[Sigma]"}], 
RowBox[{
SqrtBox["3"], " ", "m", " ", "n", " ", "sd\[Sigma]"}], 
RowBox[{
SqrtBox["3"], " ", "l", " ", "n", " ", "sd\[Sigma]"}], 
RowBox[{
FractionBox["1", "2"], " ", 
SqrtBox["3"], " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "sd\[Sigma]"}], 
RowBox[{
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}], " ", "sd\[Sigma]"}]},
{
RowBox[{
RowBox[{"-", "l"}], " ", "sp\[Sigma]R"}], 
RowBox[{
RowBox[{
RowBox[{"(", 
RowBox[{"1", "-", 
SuperscriptBox["l", "2"]}], ")"}], " ", "pp\[Pi]"}], "+", 
RowBox[{
SuperscriptBox["l", "2"], " ", "pp\[Sigma]"}]}], 
RowBox[{"l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "pp\[Pi]"}], "+", "pp\[Sigma]"}], ")"}]}], 
RowBox[{"l", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "pp\[Pi]"}], "+", "pp\[Sigma]"}], ")"}]}], 
RowBox[{
RowBox[{
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
SuperscriptBox["l", "2"]}]}], ")"}], " ", "m", " ", "pd\[Pi]"}], "+", 
RowBox[{
SqrtBox["3"], " ", 
SuperscriptBox["l", "2"], " ", "m", " ", "pd\[Sigma]"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", "2"}], " ", "l", " ", "m", " ", "n", " ", "pd\[Pi]"}], "+", 
RowBox[{
SqrtBox["3"], " ", "l", " ", "m", " ", "n", " ", "pd\[Sigma]"}]}], 
RowBox[{
RowBox[{
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
SuperscriptBox["l", "2"]}]}], ")"}], " ", "n", " ", "pd\[Pi]"}], "+", 
RowBox[{
SqrtBox["3"], " ", 
SuperscriptBox["l", "2"], " ", "n", " ", "pd\[Sigma]"}]}], 
RowBox[{
RowBox[{"l", " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"]}], ")"}], " ", "pd\[Pi]"}], "+", 
RowBox[{
FractionBox["1", "2"], " ", 
SqrtBox["3"], " ", "l", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "pd\[Sigma]"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", 
SqrtBox["3"]}], " ", "l", " ", 
SuperscriptBox["n", "2"], " ", "pd\[Pi]"}], "+", 
RowBox[{"l", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}], " ", "pd\[Sigma]"}]}]},
{
RowBox[{
RowBox[{"-", "m"}], " ", "sp\[Sigma]R"}], 
RowBox[{"l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "pp\[Pi]R"}], "+", "pp\[Sigma]R"}], ")"}]}], 
RowBox[{
RowBox[{
RowBox[{"(", 
RowBox[{"1", "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "pp\[Pi]"}], "+", 
RowBox[{
SuperscriptBox["m", "2"], " ", "pp\[Sigma]"}]}], 
RowBox[{"m", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "pp\[Pi]"}], "+", "pp\[Sigma]"}], ")"}]}], 
RowBox[{
RowBox[{"l", " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
SuperscriptBox["m", "2"]}]}], ")"}], " ", "pd\[Pi]"}], "+", 
RowBox[{
SqrtBox["3"], " ", "l", " ", 
SuperscriptBox["m", "2"], " ", "pd\[Sigma]"}]}], 
RowBox[{
RowBox[{
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
SuperscriptBox["m", "2"]}]}], ")"}], " ", "n", " ", "pd\[Pi]"}], "+", 
RowBox[{
SqrtBox["3"], " ", 
SuperscriptBox["m", "2"], " ", "n", " ", "pd\[Sigma]"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", "2"}], " ", "l", " ", "m", " ", "n", " ", "pd\[Pi]"}], "+", 
RowBox[{
SqrtBox["3"], " ", "l", " ", "m", " ", "n", " ", "pd\[Sigma]"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", "m"}], " ", 
RowBox[{"(", 
RowBox[{"1", "+", 
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "pd\[Pi]"}], "+", 
RowBox[{
FractionBox["1", "2"], " ", 
SqrtBox["3"], " ", "m", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "pd\[Sigma]"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", 
SqrtBox["3"]}], " ", "m", " ", 
SuperscriptBox["n", "2"], " ", "pd\[Pi]"}], "+", 
RowBox[{"m", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}], " ", "pd\[Sigma]"}]}]},
{
RowBox[{
RowBox[{"-", "n"}], " ", "sp\[Sigma]R"}], 
RowBox[{"l", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "pp\[Pi]R"}], "+", "pp\[Sigma]R"}], ")"}]}], 
RowBox[{"m", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "pp\[Pi]R"}], "+", "pp\[Sigma]R"}], ")"}]}], 
RowBox[{
RowBox[{
RowBox[{"(", 
RowBox[{"1", "-", 
SuperscriptBox["n", "2"]}], ")"}], " ", "pp\[Pi]"}], "+", 
RowBox[{
SuperscriptBox["n", "2"], " ", "pp\[Sigma]"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", "2"}], " ", "l", " ", "m", " ", "n", " ", "pd\[Pi]"}], "+", 
RowBox[{
SqrtBox["3"], " ", "l", " ", "m", " ", "n", " ", "pd\[Sigma]"}]}], 
RowBox[{
RowBox[{"m", " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
SuperscriptBox["n", "2"]}]}], ")"}], " ", "pd\[Pi]"}], "+", 
RowBox[{
SqrtBox["3"], " ", "m", " ", 
SuperscriptBox["n", "2"], " ", "pd\[Sigma]"}]}], 
RowBox[{
RowBox[{"l", " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
SuperscriptBox["n", "2"]}]}], ")"}], " ", "pd\[Pi]"}], "+", 
RowBox[{
SqrtBox["3"], " ", "l", " ", 
SuperscriptBox["n", "2"], " ", "pd\[Sigma]"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], " ", "n", " ", "pd\[Pi]"}], "+", 
RowBox[{
FractionBox["1", "2"], " ", 
SqrtBox["3"], " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n", " ", "pd\[Sigma]"}]}], 
RowBox[{
RowBox[{
SqrtBox["3"], " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n", " ", "pd\[Pi]"}], "+", 
RowBox[{"n", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["n", "2"], "-", 
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"]}], ")"}]}]}], ")"}], " ", "pd\[Sigma]"}]}]},
{
RowBox[{
SqrtBox["3"], " ", "l", " ", "m", " ", "sd\[Sigma]R"}], 
RowBox[{
RowBox[{
RowBox[{"-", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
SuperscriptBox["l", "2"]}]}], ")"}]}], " ", "m", " ", "pd\[Pi]R"}], "-", 
RowBox[{
SqrtBox["3"], " ", 
SuperscriptBox["l", "2"], " ", "m", " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", "l"}], " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
SuperscriptBox["m", "2"]}]}], ")"}], " ", "pd\[Pi]R"}], "-", 
RowBox[{
SqrtBox["3"], " ", "l", " ", 
SuperscriptBox["m", "2"], " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{"2", " ", "l", " ", "m", " ", "n", " ", "pd\[Pi]R"}], "-", 
RowBox[{
SqrtBox["3"], " ", "l", " ", "m", " ", "n", " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{"3", " ", "dd\[Sigma]", " ", 
SuperscriptBox["l", "2"], " ", 
SuperscriptBox["m", "2"]}], "+", 
RowBox[{"dd\[Pi]", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"], "-", 
RowBox[{"4", " ", 
SuperscriptBox["l", "2"], " ", 
SuperscriptBox["m", "2"]}]}], ")"}]}], "+", 
RowBox[{"dd\[Delta]", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
SuperscriptBox["l", "2"], " ", 
SuperscriptBox["m", "2"]}], "+", 
SuperscriptBox["n", "2"]}], ")"}]}]}], 
RowBox[{
RowBox[{"3", " ", "dd\[Sigma]", " ", "l", " ", 
SuperscriptBox["m", "2"], " ", "n"}], "+", 
RowBox[{"dd\[Pi]", " ", "l", " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"4", " ", 
SuperscriptBox["m", "2"]}]}], ")"}], " ", "n"}], "+", 
RowBox[{"dd\[Delta]", " ", "l", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "1"}], "+", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n"}]}], 
RowBox[{
RowBox[{"3", " ", "dd\[Sigma]", " ", 
SuperscriptBox["l", "2"], " ", "m", " ", "n"}], "+", 
RowBox[{"dd\[Pi]", " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"4", " ", 
SuperscriptBox["l", "2"]}]}], ")"}], " ", "m", " ", "n"}], "+", 
RowBox[{"dd\[Delta]", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "1"}], "+", 
SuperscriptBox["l", "2"]}], ")"}], " ", "m", " ", "n"}]}], 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", "dd\[Delta]", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
RowBox[{
FractionBox["3", "2"], " ", "dd\[Sigma]", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
RowBox[{"2", " ", "dd\[Pi]", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "+", 
SuperscriptBox["m", "2"]}], ")"}]}]}], 
RowBox[{
RowBox[{
RowBox[{"-", "2"}], " ", 
SqrtBox["3"], " ", "dd\[Pi]", " ", "l", " ", "m", " ", 
SuperscriptBox["n", "2"]}], "+", 
RowBox[{
FractionBox["1", "2"], " ", 
SqrtBox["3"], " ", "dd\[Delta]", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{"1", "+", 
SuperscriptBox["n", "2"]}], ")"}]}], "+", 
RowBox[{
SqrtBox["3"], " ", "dd\[Sigma]", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}]}]}]},
{
RowBox[{
SqrtBox["3"], " ", "m", " ", "n", " ", "sd\[Sigma]R"}], 
RowBox[{
RowBox[{"2", " ", "l", " ", "m", " ", "n", " ", "pd\[Pi]R"}], "-", 
RowBox[{
SqrtBox["3"], " ", "l", " ", "m", " ", "n", " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
SuperscriptBox["m", "2"]}]}], ")"}]}], " ", "n", " ", "pd\[Pi]R"}], "-", 
RowBox[{
SqrtBox["3"], " ", 
SuperscriptBox["m", "2"], " ", "n", " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", "m"}], " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
SuperscriptBox["n", "2"]}]}], ")"}], " ", "pd\[Pi]R"}], "-", 
RowBox[{
SqrtBox["3"], " ", "m", " ", 
SuperscriptBox["n", "2"], " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{"3", " ", "dd\[Sigma]R", " ", "l", " ", 
SuperscriptBox["m", "2"], " ", "n"}], "+", 
RowBox[{"dd\[Pi]R", " ", "l", " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"4", " ", 
SuperscriptBox["m", "2"]}]}], ")"}], " ", "n"}], "+", 
RowBox[{"dd\[Delta]", " ", "l", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "1"}], "+", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n"}]}], 
RowBox[{
RowBox[{"3", " ", "dd\[Sigma]", " ", 
SuperscriptBox["m", "2"], " ", 
SuperscriptBox["n", "2"]}], "+", 
RowBox[{"dd\[Pi]", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["m", "2"], "+", 
SuperscriptBox["n", "2"], "-", 
RowBox[{"4", " ", 
SuperscriptBox["m", "2"], " ", 
SuperscriptBox["n", "2"]}]}], ")"}]}], "+", 
RowBox[{"dd\[Delta]", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
RowBox[{
SuperscriptBox["m", "2"], " ", 
SuperscriptBox["n", "2"]}]}], ")"}]}]}], 
RowBox[{
RowBox[{"3", " ", "dd\[Sigma]", " ", "l", " ", "m", " ", 
SuperscriptBox["n", "2"]}], "+", 
RowBox[{"dd\[Pi]", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"4", " ", 
SuperscriptBox["n", "2"]}]}], ")"}]}], "+", 
RowBox[{"dd\[Delta]", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "1"}], "+", 
SuperscriptBox["n", "2"]}], ")"}]}]}], 
RowBox[{
RowBox[{
FractionBox["3", "2"], " ", "dd\[Sigma]", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n"}], "+", 
RowBox[{"dd\[Delta]", " ", "m", " ", 
RowBox[{"(", 
RowBox[{"1", "+", 
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}]}]}], ")"}], " ", "n"}], "-", 
RowBox[{"dd\[Pi]", " ", "m", " ", 
RowBox[{"(", 
RowBox[{"1", "+", 
RowBox[{"2", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}]}]}], ")"}], " ", "n"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", 
FractionBox["1", "2"]}], " ", 
SqrtBox["3"], " ", "dd\[Delta]", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n"}], "+", 
RowBox[{
SqrtBox["3"], " ", "dd\[Pi]", " ", "m", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"], "-", 
SuperscriptBox["n", "2"]}], ")"}]}], "+", 
RowBox[{
SqrtBox["3"], " ", "dd\[Sigma]", "  ", "m", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}]}]}]},
{
RowBox[{
SqrtBox["3"], " ", "l", " ", "n", " ", "sd\[Sigma]R"}], 
RowBox[{
RowBox[{
RowBox[{"-", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
SuperscriptBox["l", "2"]}]}], ")"}]}], " ", "n", " ", "pd\[Pi]R"}], "-", 
RowBox[{
SqrtBox["3"], " ", 
SuperscriptBox["l", "2"], " ", "n", " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{"2", " ", "l", " ", "m", " ", "n", " ", "pd\[Pi]R"}], "-", 
RowBox[{
SqrtBox["3"], " ", "l", " ", "m", " ", "n", " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", "l"}], " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
SuperscriptBox["n", "2"]}]}], ")"}], " ", "pd\[Pi]R"}], "-", 
RowBox[{
SqrtBox["3"], " ", "l", " ", 
SuperscriptBox["n", "2"], " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{"3", " ", "dd\[Sigma]R", " ", 
SuperscriptBox["l", "2"], " ", "m", " ", "n"}], "+", 
RowBox[{"dd\[Pi]R", " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"4", " ", 
SuperscriptBox["l", "2"]}]}], ")"}], " ", "m", " ", "n"}], "+", 
RowBox[{"dd\[Delta]", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "1"}], "+", 
SuperscriptBox["l", "2"]}], ")"}], " ", "m", " ", "n"}]}], 
RowBox[{
RowBox[{"3", " ", "dd\[Sigma]R", " ", "l", " ", "m", " ", 
SuperscriptBox["n", "2"]}], "+", 
RowBox[{"dd\[Pi]R", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"4", " ", 
SuperscriptBox["n", "2"]}]}], ")"}]}], "+", 
RowBox[{"dd\[Delta]", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "1"}], "+", 
SuperscriptBox["n", "2"]}], ")"}]}]}], 
RowBox[{
RowBox[{"3", " ", "dd\[Sigma]", " ", 
SuperscriptBox["l", "2"], " ", 
SuperscriptBox["n", "2"]}], "+", 
RowBox[{"dd\[Pi]", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["n", "2"], "-", 
RowBox[{"4", " ", 
SuperscriptBox["l", "2"], " ", 
SuperscriptBox["n", "2"]}]}], ")"}]}], "+", 
RowBox[{"dd\[Delta]", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["m", "2"], "+", 
RowBox[{
SuperscriptBox["l", "2"], " ", 
SuperscriptBox["n", "2"]}]}], ")"}]}]}], 
RowBox[{
RowBox[{
FractionBox["3", "2"], " ", "dd\[Sigma]", " ", "l", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n"}], "+", 
RowBox[{"dd\[Pi]", " ", "l", " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}]}]}], ")"}], " ", "n"}], "-", 
RowBox[{"dd\[Delta]", " ", "l", " ", 
RowBox[{"(", 
RowBox[{"1", "+", 
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "+", 
SuperscriptBox["m", "2"]}], ")"}]}]}], ")"}], " ", "n"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", 
FractionBox["1", "2"]}], " ", 
SqrtBox["3"], " ", "dd\[Delta]", " ", "l", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n"}], "+", 
RowBox[{
SqrtBox["3"], " ", "dd\[Pi]", " ", "l", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"], "-", 
SuperscriptBox["n", "2"]}], ")"}]}], "+", 
RowBox[{
SqrtBox["3"], " ", "dd\[Sigma]", " ", "l", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}]}]}]},
{
RowBox[{
FractionBox["1", "2"], " ", 
SqrtBox["3"], " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "sd\[Sigma]R"}], 
RowBox[{
RowBox[{
RowBox[{"-", "l"}], " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"]}], ")"}], " ", "pd\[Pi]R"}], "-", 
RowBox[{
FractionBox["1", "2"], " ", 
SqrtBox["3"], " ", "l", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{"m", " ", 
RowBox[{"(", 
RowBox[{"1", "+", 
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "pd\[Pi]R"}], "-", 
RowBox[{
FractionBox["1", "2"], " ", 
SqrtBox["3"], " ", "m", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n", " ", "pd\[Pi]R"}], "-", 
RowBox[{
FractionBox["1", "2"], " ", 
SqrtBox["3"], " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n", " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", "dd\[Delta]", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
RowBox[{
FractionBox["3", "2"], " ", "dd\[Sigma]R", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
RowBox[{"2", " ", "dd\[Pi]R", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "+", 
SuperscriptBox["m", "2"]}], ")"}]}]}], 
RowBox[{
RowBox[{
FractionBox["3", "2"], " ", "dd\[Sigma]R", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n"}], "+", 
RowBox[{"dd\[Delta]", " ", "m", " ", 
RowBox[{"(", 
RowBox[{"1", "+", 
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}]}]}], ")"}], " ", "n"}], "-", 
RowBox[{"dd\[Pi]R", " ", "m", " ", 
RowBox[{"(", 
RowBox[{"1", "+", 
RowBox[{"2", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}]}]}], ")"}], " ", "n"}]}], 
RowBox[{
RowBox[{
FractionBox["3", "2"], " ", "dd\[Sigma]R", " ", "l", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n"}], "+", 
RowBox[{"dd\[Pi]R", " ", "l", " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}]}]}], ")"}], " ", "n"}], "-", 
RowBox[{"dd\[Delta]", " ", "l", " ", 
RowBox[{"(", 
RowBox[{"1", "+", 
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "+", 
SuperscriptBox["m", "2"]}], ")"}]}]}], ")"}], " ", "n"}]}], 
RowBox[{
RowBox[{
FractionBox["3", "4"], " ", "dd\[Sigma]", " ", 
SuperscriptBox[
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], "2"]}], "+", 
RowBox[{"dd\[Pi]", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"], "-", 
SuperscriptBox[
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], "2"]}], ")"}]}], "+", 
RowBox[{"dd\[Delta]", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "4"], " ", 
SuperscriptBox[
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], "2"]}], "+", 
SuperscriptBox["n", "2"]}], ")"}]}]}], 
RowBox[{
RowBox[{
SqrtBox["3"], " ", "dd\[Pi]", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "+", 
SuperscriptBox["m", "2"]}], ")"}], " ", 
SuperscriptBox["n", "2"]}], "+", 
RowBox[{
FractionBox["1", "4"], " ", 
SqrtBox["3"], " ", "dd\[Delta]", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", 
RowBox[{"(", 
RowBox[{"1", "+", 
SuperscriptBox["n", "2"]}], ")"}]}], "+", 
RowBox[{
FractionBox["1", "2"], " ", 
SqrtBox["3"], " ", "dd\[Sigma]", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}]}]}]},
{
RowBox[{
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}], " ", "sd\[Sigma]R"}], 
RowBox[{
RowBox[{
SqrtBox["3"], " ", "l", " ", 
SuperscriptBox["n", "2"], " ", "pd\[Pi]R"}], "-", 
RowBox[{"l", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}], " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{
SqrtBox["3"], " ", "m", " ", 
SuperscriptBox["n", "2"], " ", "pd\[Pi]R"}], "-", 
RowBox[{"m", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}], " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", 
SqrtBox["3"]}], " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n", " ", "pd\[Pi]R"}], "-", 
RowBox[{"n", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}], " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", "2"}], " ", 
SqrtBox["3"], " ", "dd\[Pi]R", " ", "l", " ", "m", " ", 
SuperscriptBox["n", "2"]}], "+", 
RowBox[{
FractionBox["1", "2"], " ", 
SqrtBox["3"], " ", "dd\[Delta]", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{"1", "+", 
SuperscriptBox["n", "2"]}], ")"}]}], "+", 
RowBox[{
SqrtBox["3"], " ", "dd\[Sigma]R", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}]}]}], 
RowBox[{
RowBox[{
RowBox[{"-", 
FractionBox["1", "2"]}], " ", 
SqrtBox["3"], " ", "dd\[Delta]", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n"}], "+", 
RowBox[{
SqrtBox["3"], " ", "dd\[Pi]R", " ", "m", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"], "-", 
SuperscriptBox["n", "2"]}], ")"}]}], "+", 
RowBox[{
SqrtBox["3"], " ", "dd\[Sigma]R", " ", "m", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}]}]}], 
RowBox[{
RowBox[{
RowBox[{"-", 
FractionBox["1", "2"]}], " ", 
SqrtBox["3"], " ", "dd\[Delta]", " ", "l", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n"}], "+", 
RowBox[{
SqrtBox["3"], " ", "dd\[Pi]R", " ", "l", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"], "-", 
SuperscriptBox["n", "2"]}], ")"}]}], "+", 
RowBox[{
SqrtBox["3"], " ", "dd\[Sigma]R", " ", "l", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}]}]}], 
RowBox[{
RowBox[{
SqrtBox["3"], " ", "dd\[Pi]R", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "+", 
SuperscriptBox["m", "2"]}], ")"}], " ", 
SuperscriptBox["n", "2"]}], "+", 
RowBox[{
FractionBox["1", "4"], " ", 
SqrtBox["3"], " ", "dd\[Delta]", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", 
RowBox[{"(", 
RowBox[{"1", "+", 
SuperscriptBox["n", "2"]}], ")"}]}], "+", 
RowBox[{
FractionBox["1", "2"], " ", 
SqrtBox["3"], " ", "dd\[Sigma]R", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}]}]}], 
RowBox[{
RowBox[{
SuperscriptBox[
RowBox[{"(", 
RowBox[{
SuperscriptBox["n", "2"], "-", 
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"]}], ")"}]}]}], ")"}], "2"], "dd\[Sigma]"}], " ", "+", 
RowBox[{"3", " ", 
SuperscriptBox["n", "2"], 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"]}], ")"}], " ", "dd\[Pi]"}], " ", "+", 
RowBox[{
FractionBox["3", "4"], "  ", 
SuperscriptBox[
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"]}], ")"}], "2"], "dd\[Delta]"}]}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\);
];
];
];

(*all real, transpose is ok*)
(*pureHint=RU+RU\[Transpose]-RU IdentityMatrix[9]//N;*)
(*t2=(ret=pureHint[[1;;MAXORB[a1],1;;MAXORB[a2]]]//N;)//AbsoluteTiming;*)
pureHint[[1;;MAXORB[a1],1;;MAXORB[a2]]]//N
(*IdentityMatrix[9][[1;;MAXORB[a1],1;;MAXORB[a2]]]//N*)
]

SlaterOver[R_,a1_,a2_]:=Module[{bond=a1~~"-"~~a2,
bondR=a2~~"-"~~a1,
r=Normalize[R],
rr=Norm[R],
maxorb2=Max[MAXORB/@{a1,a2}],
l,m,n,ss\[Sigma],sp\[Sigma],sd\[Sigma],pp\[Sigma],pd\[Sigma],dd\[Sigma],pp\[Pi],pd\[Pi],dd\[Pi],dd\[Delta],ss\[Sigma]R,sp\[Sigma]R,sd\[Sigma]R,pp\[Sigma]R,pd\[Sigma]R,dd\[Sigma]R,pp\[Pi]R,pd\[Pi]R,dd\[Pi]R,dd\[Delta]R,pureHint},

{l,m,n}=r;

ss\[Sigma]=S[bond~~".ss\[Sigma]"][rr];
ss\[Sigma]R=S[bondR~~".ss\[Sigma]"][rr];
If[maxorb2>1,
sp\[Sigma]=S[bond~~".sp\[Sigma]"][rr];
pp\[Sigma]=S[bond~~".pp\[Sigma]"][rr];
pp\[Pi]=S[bond~~".pp\[Pi]"][rr];

sp\[Sigma]R=S[bondR~~".sp\[Sigma]"][rr];
pp\[Sigma]R=S[bondR~~".pp\[Sigma]"][rr];
pp\[Pi]R=S[bondR~~".pp\[Pi]"][rr];
];
If[maxorb2>4,
sd\[Sigma]=S[bond~~".sd\[Sigma]"][rr];
dd\[Sigma]=S[bond~~".dd\[Sigma]"][rr];
pd\[Sigma]=S[bond~~".pd\[Sigma]"][rr];
pd\[Pi]=S[bond~~".pd\[Pi]"][rr];
dd\[Pi]=S[bond~~".dd\[Pi]"][rr];
dd\[Delta]=S[bond~~".dd\[Delta]"][rr];
sd\[Sigma]R=S[bondR~~".sd\[Sigma]"][rr];
pd\[Sigma]R=S[bondR~~".pd\[Sigma]"][rr];
dd\[Sigma]R=S[bondR~~".dd\[Sigma]"][rr];
pd\[Pi]R=S[bondR~~".pd\[Pi]"][rr];
dd\[Pi]R=S[bondR~~".dd\[Pi]"][rr];
dd\[Delta]R=S[bondR~~".dd\[Delta]"][rr];
];

(*s x y z xy yz zx x^2-y^2 3z^2-r^2*)
If[maxorb2==1,pureHint=pureHint={{ss\[Sigma]}},
If[maxorb2==4,pureHint=\!\(\*
TagBox[
TagBox[
TagBox[
RowBox[{"(", GridBox[{
{"ss\[Sigma]", 
RowBox[{"l", " ", "sp\[Sigma]"}], 
RowBox[{"m", " ", "sp\[Sigma]"}], 
RowBox[{"n", " ", "sp\[Sigma]"}]},
{
RowBox[{
RowBox[{"-", "l"}], " ", "sp\[Sigma]R"}], 
RowBox[{
RowBox[{
RowBox[{"(", 
RowBox[{"1", "-", 
SuperscriptBox["l", "2"]}], ")"}], " ", "pp\[Pi]"}], "+", 
RowBox[{
SuperscriptBox["l", "2"], " ", "pp\[Sigma]"}]}], 
RowBox[{"l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "pp\[Pi]"}], "+", "pp\[Sigma]"}], ")"}]}], 
RowBox[{"l", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "pp\[Pi]"}], "+", "pp\[Sigma]"}], ")"}]}]},
{
RowBox[{
RowBox[{"-", "m"}], " ", "sp\[Sigma]R"}], 
RowBox[{"l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "pp\[Pi]R"}], "+", "pp\[Sigma]R"}], ")"}]}], 
RowBox[{
RowBox[{
RowBox[{"(", 
RowBox[{"1", "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "pp\[Pi]"}], "+", 
RowBox[{
SuperscriptBox["m", "2"], " ", "pp\[Sigma]"}]}], 
RowBox[{"m", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "pp\[Pi]"}], "+", "pp\[Sigma]"}], ")"}]}]},
{
RowBox[{
RowBox[{"-", "n"}], " ", "sp\[Sigma]R"}], 
RowBox[{"l", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "pp\[Pi]R"}], "+", "pp\[Sigma]R"}], ")"}]}], 
RowBox[{"m", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "pp\[Pi]R"}], "+", "pp\[Sigma]R"}], ")"}]}], 
RowBox[{
RowBox[{
RowBox[{"(", 
RowBox[{"1", "-", 
SuperscriptBox["n", "2"]}], ")"}], " ", "pp\[Pi]"}], "+", 
RowBox[{
SuperscriptBox["n", "2"], " ", "pp\[Sigma]"}]}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\),
If[maxorb2==9,pureHint=\!\(\*
TagBox[
TagBox[
TagBox[
RowBox[{"(", GridBox[{
{"ss\[Sigma]", 
RowBox[{"l", " ", "sp\[Sigma]"}], 
RowBox[{"m", " ", "sp\[Sigma]"}], 
RowBox[{"n", " ", "sp\[Sigma]"}], 
RowBox[{
SqrtBox["3"], " ", "l", " ", "m", " ", "sd\[Sigma]"}], 
RowBox[{
SqrtBox["3"], " ", "m", " ", "n", " ", "sd\[Sigma]"}], 
RowBox[{
SqrtBox["3"], " ", "l", " ", "n", " ", "sd\[Sigma]"}], 
RowBox[{
FractionBox["1", "2"], " ", 
SqrtBox["3"], " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "sd\[Sigma]"}], 
RowBox[{
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}], " ", "sd\[Sigma]"}]},
{
RowBox[{
RowBox[{"-", "l"}], " ", "sp\[Sigma]R"}], 
RowBox[{
RowBox[{
RowBox[{"(", 
RowBox[{"1", "-", 
SuperscriptBox["l", "2"]}], ")"}], " ", "pp\[Pi]"}], "+", 
RowBox[{
SuperscriptBox["l", "2"], " ", "pp\[Sigma]"}]}], 
RowBox[{"l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "pp\[Pi]"}], "+", "pp\[Sigma]"}], ")"}]}], 
RowBox[{"l", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "pp\[Pi]"}], "+", "pp\[Sigma]"}], ")"}]}], 
RowBox[{
RowBox[{
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
SuperscriptBox["l", "2"]}]}], ")"}], " ", "m", " ", "pd\[Pi]"}], "+", 
RowBox[{
SqrtBox["3"], " ", 
SuperscriptBox["l", "2"], " ", "m", " ", "pd\[Sigma]"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", "2"}], " ", "l", " ", "m", " ", "n", " ", "pd\[Pi]"}], "+", 
RowBox[{
SqrtBox["3"], " ", "l", " ", "m", " ", "n", " ", "pd\[Sigma]"}]}], 
RowBox[{
RowBox[{
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
SuperscriptBox["l", "2"]}]}], ")"}], " ", "n", " ", "pd\[Pi]"}], "+", 
RowBox[{
SqrtBox["3"], " ", 
SuperscriptBox["l", "2"], " ", "n", " ", "pd\[Sigma]"}]}], 
RowBox[{
RowBox[{"l", " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"]}], ")"}], " ", "pd\[Pi]"}], "+", 
RowBox[{
FractionBox["1", "2"], " ", 
SqrtBox["3"], " ", "l", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "pd\[Sigma]"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", 
SqrtBox["3"]}], " ", "l", " ", 
SuperscriptBox["n", "2"], " ", "pd\[Pi]"}], "+", 
RowBox[{"l", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}], " ", "pd\[Sigma]"}]}]},
{
RowBox[{
RowBox[{"-", "m"}], " ", "sp\[Sigma]R"}], 
RowBox[{"l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "pp\[Pi]R"}], "+", "pp\[Sigma]R"}], ")"}]}], 
RowBox[{
RowBox[{
RowBox[{"(", 
RowBox[{"1", "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "pp\[Pi]"}], "+", 
RowBox[{
SuperscriptBox["m", "2"], " ", "pp\[Sigma]"}]}], 
RowBox[{"m", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "pp\[Pi]"}], "+", "pp\[Sigma]"}], ")"}]}], 
RowBox[{
RowBox[{"l", " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
SuperscriptBox["m", "2"]}]}], ")"}], " ", "pd\[Pi]"}], "+", 
RowBox[{
SqrtBox["3"], " ", "l", " ", 
SuperscriptBox["m", "2"], " ", "pd\[Sigma]"}]}], 
RowBox[{
RowBox[{
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
SuperscriptBox["m", "2"]}]}], ")"}], " ", "n", " ", "pd\[Pi]"}], "+", 
RowBox[{
SqrtBox["3"], " ", 
SuperscriptBox["m", "2"], " ", "n", " ", "pd\[Sigma]"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", "2"}], " ", "l", " ", "m", " ", "n", " ", "pd\[Pi]"}], "+", 
RowBox[{
SqrtBox["3"], " ", "l", " ", "m", " ", "n", " ", "pd\[Sigma]"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", "m"}], " ", 
RowBox[{"(", 
RowBox[{"1", "+", 
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "pd\[Pi]"}], "+", 
RowBox[{
FractionBox["1", "2"], " ", 
SqrtBox["3"], " ", "m", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "pd\[Sigma]"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", 
SqrtBox["3"]}], " ", "m", " ", 
SuperscriptBox["n", "2"], " ", "pd\[Pi]"}], "+", 
RowBox[{"m", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}], " ", "pd\[Sigma]"}]}]},
{
RowBox[{
RowBox[{"-", "n"}], " ", "sp\[Sigma]R"}], 
RowBox[{"l", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "pp\[Pi]R"}], "+", "pp\[Sigma]R"}], ")"}]}], 
RowBox[{"m", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "pp\[Pi]R"}], "+", "pp\[Sigma]R"}], ")"}]}], 
RowBox[{
RowBox[{
RowBox[{"(", 
RowBox[{"1", "-", 
SuperscriptBox["n", "2"]}], ")"}], " ", "pp\[Pi]"}], "+", 
RowBox[{
SuperscriptBox["n", "2"], " ", "pp\[Sigma]"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", "2"}], " ", "l", " ", "m", " ", "n", " ", "pd\[Pi]"}], "+", 
RowBox[{
SqrtBox["3"], " ", "l", " ", "m", " ", "n", " ", "pd\[Sigma]"}]}], 
RowBox[{
RowBox[{"m", " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
SuperscriptBox["n", "2"]}]}], ")"}], " ", "pd\[Pi]"}], "+", 
RowBox[{
SqrtBox["3"], " ", "m", " ", 
SuperscriptBox["n", "2"], " ", "pd\[Sigma]"}]}], 
RowBox[{
RowBox[{"l", " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
SuperscriptBox["n", "2"]}]}], ")"}], " ", "pd\[Pi]"}], "+", 
RowBox[{
SqrtBox["3"], " ", "l", " ", 
SuperscriptBox["n", "2"], " ", "pd\[Sigma]"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], " ", "n", " ", "pd\[Pi]"}], "+", 
RowBox[{
FractionBox["1", "2"], " ", 
SqrtBox["3"], " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n", " ", "pd\[Sigma]"}]}], 
RowBox[{
RowBox[{
SqrtBox["3"], " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n", " ", "pd\[Pi]"}], "+", 
RowBox[{"n", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["n", "2"], "-", 
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"]}], ")"}]}]}], ")"}], " ", "pd\[Sigma]"}]}]},
{
RowBox[{
SqrtBox["3"], " ", "l", " ", "m", " ", "sd\[Sigma]R"}], 
RowBox[{
RowBox[{
RowBox[{"-", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
SuperscriptBox["l", "2"]}]}], ")"}]}], " ", "m", " ", "pd\[Pi]R"}], "-", 
RowBox[{
SqrtBox["3"], " ", 
SuperscriptBox["l", "2"], " ", "m", " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", "l"}], " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
SuperscriptBox["m", "2"]}]}], ")"}], " ", "pd\[Pi]R"}], "-", 
RowBox[{
SqrtBox["3"], " ", "l", " ", 
SuperscriptBox["m", "2"], " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{"2", " ", "l", " ", "m", " ", "n", " ", "pd\[Pi]R"}], "-", 
RowBox[{
SqrtBox["3"], " ", "l", " ", "m", " ", "n", " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{"3", " ", "dd\[Sigma]", " ", 
SuperscriptBox["l", "2"], " ", 
SuperscriptBox["m", "2"]}], "+", 
RowBox[{"dd\[Pi]", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"], "-", 
RowBox[{"4", " ", 
SuperscriptBox["l", "2"], " ", 
SuperscriptBox["m", "2"]}]}], ")"}]}], "+", 
RowBox[{"dd\[Delta]", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
SuperscriptBox["l", "2"], " ", 
SuperscriptBox["m", "2"]}], "+", 
SuperscriptBox["n", "2"]}], ")"}]}]}], 
RowBox[{
RowBox[{"3", " ", "dd\[Sigma]", " ", "l", " ", 
SuperscriptBox["m", "2"], " ", "n"}], "+", 
RowBox[{"dd\[Pi]", " ", "l", " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"4", " ", 
SuperscriptBox["m", "2"]}]}], ")"}], " ", "n"}], "+", 
RowBox[{"dd\[Delta]", " ", "l", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "1"}], "+", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n"}]}], 
RowBox[{
RowBox[{"3", " ", "dd\[Sigma]", " ", 
SuperscriptBox["l", "2"], " ", "m", " ", "n"}], "+", 
RowBox[{"dd\[Pi]", " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"4", " ", 
SuperscriptBox["l", "2"]}]}], ")"}], " ", "m", " ", "n"}], "+", 
RowBox[{"dd\[Delta]", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "1"}], "+", 
SuperscriptBox["l", "2"]}], ")"}], " ", "m", " ", "n"}]}], 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", "dd\[Delta]", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
RowBox[{
FractionBox["3", "2"], " ", "dd\[Sigma]", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
RowBox[{"2", " ", "dd\[Pi]", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "+", 
SuperscriptBox["m", "2"]}], ")"}]}]}], 
RowBox[{
RowBox[{
RowBox[{"-", "2"}], " ", 
SqrtBox["3"], " ", "dd\[Pi]", " ", "l", " ", "m", " ", 
SuperscriptBox["n", "2"]}], "+", 
RowBox[{
FractionBox["1", "2"], " ", 
SqrtBox["3"], " ", "dd\[Delta]", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{"1", "+", 
SuperscriptBox["n", "2"]}], ")"}]}], "+", 
RowBox[{
SqrtBox["3"], " ", "dd\[Sigma]", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}]}]}]},
{
RowBox[{
SqrtBox["3"], " ", "m", " ", "n", " ", "sd\[Sigma]R"}], 
RowBox[{
RowBox[{"2", " ", "l", " ", "m", " ", "n", " ", "pd\[Pi]R"}], "-", 
RowBox[{
SqrtBox["3"], " ", "l", " ", "m", " ", "n", " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
SuperscriptBox["m", "2"]}]}], ")"}]}], " ", "n", " ", "pd\[Pi]R"}], "-", 
RowBox[{
SqrtBox["3"], " ", 
SuperscriptBox["m", "2"], " ", "n", " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", "m"}], " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
SuperscriptBox["n", "2"]}]}], ")"}], " ", "pd\[Pi]R"}], "-", 
RowBox[{
SqrtBox["3"], " ", "m", " ", 
SuperscriptBox["n", "2"], " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{"3", " ", "dd\[Sigma]R", " ", "l", " ", 
SuperscriptBox["m", "2"], " ", "n"}], "+", 
RowBox[{"dd\[Pi]R", " ", "l", " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"4", " ", 
SuperscriptBox["m", "2"]}]}], ")"}], " ", "n"}], "+", 
RowBox[{"dd\[Delta]", " ", "l", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "1"}], "+", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n"}]}], 
RowBox[{
RowBox[{"3", " ", "dd\[Sigma]", " ", 
SuperscriptBox["m", "2"], " ", 
SuperscriptBox["n", "2"]}], "+", 
RowBox[{"dd\[Pi]", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["m", "2"], "+", 
SuperscriptBox["n", "2"], "-", 
RowBox[{"4", " ", 
SuperscriptBox["m", "2"], " ", 
SuperscriptBox["n", "2"]}]}], ")"}]}], "+", 
RowBox[{"dd\[Delta]", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
RowBox[{
SuperscriptBox["m", "2"], " ", 
SuperscriptBox["n", "2"]}]}], ")"}]}]}], 
RowBox[{
RowBox[{"3", " ", "dd\[Sigma]", " ", "l", " ", "m", " ", 
SuperscriptBox["n", "2"]}], "+", 
RowBox[{"dd\[Pi]", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"4", " ", 
SuperscriptBox["n", "2"]}]}], ")"}]}], "+", 
RowBox[{"dd\[Delta]", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "1"}], "+", 
SuperscriptBox["n", "2"]}], ")"}]}]}], 
RowBox[{
RowBox[{
FractionBox["3", "2"], " ", "dd\[Sigma]", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n"}], "+", 
RowBox[{"dd\[Delta]", " ", "m", " ", 
RowBox[{"(", 
RowBox[{"1", "+", 
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}]}]}], ")"}], " ", "n"}], "-", 
RowBox[{"dd\[Pi]", " ", "m", " ", 
RowBox[{"(", 
RowBox[{"1", "+", 
RowBox[{"2", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}]}]}], ")"}], " ", "n"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", 
FractionBox["1", "2"]}], " ", 
SqrtBox["3"], " ", "dd\[Delta]", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n"}], "+", 
RowBox[{
SqrtBox["3"], " ", "dd\[Pi]", " ", "m", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"], "-", 
SuperscriptBox["n", "2"]}], ")"}]}], "+", 
RowBox[{
SqrtBox["3"], " ", "dd\[Sigma]", "  ", "m", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}]}]}]},
{
RowBox[{
SqrtBox["3"], " ", "l", " ", "n", " ", "sd\[Sigma]R"}], 
RowBox[{
RowBox[{
RowBox[{"-", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
SuperscriptBox["l", "2"]}]}], ")"}]}], " ", "n", " ", "pd\[Pi]R"}], "-", 
RowBox[{
SqrtBox["3"], " ", 
SuperscriptBox["l", "2"], " ", "n", " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{"2", " ", "l", " ", "m", " ", "n", " ", "pd\[Pi]R"}], "-", 
RowBox[{
SqrtBox["3"], " ", "l", " ", "m", " ", "n", " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", "l"}], " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
SuperscriptBox["n", "2"]}]}], ")"}], " ", "pd\[Pi]R"}], "-", 
RowBox[{
SqrtBox["3"], " ", "l", " ", 
SuperscriptBox["n", "2"], " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{"3", " ", "dd\[Sigma]R", " ", 
SuperscriptBox["l", "2"], " ", "m", " ", "n"}], "+", 
RowBox[{"dd\[Pi]R", " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"4", " ", 
SuperscriptBox["l", "2"]}]}], ")"}], " ", "m", " ", "n"}], "+", 
RowBox[{"dd\[Delta]", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "1"}], "+", 
SuperscriptBox["l", "2"]}], ")"}], " ", "m", " ", "n"}]}], 
RowBox[{
RowBox[{"3", " ", "dd\[Sigma]R", " ", "l", " ", "m", " ", 
SuperscriptBox["n", "2"]}], "+", 
RowBox[{"dd\[Pi]R", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"4", " ", 
SuperscriptBox["n", "2"]}]}], ")"}]}], "+", 
RowBox[{"dd\[Delta]", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", "1"}], "+", 
SuperscriptBox["n", "2"]}], ")"}]}]}], 
RowBox[{
RowBox[{"3", " ", "dd\[Sigma]", " ", 
SuperscriptBox["l", "2"], " ", 
SuperscriptBox["n", "2"]}], "+", 
RowBox[{"dd\[Pi]", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["n", "2"], "-", 
RowBox[{"4", " ", 
SuperscriptBox["l", "2"], " ", 
SuperscriptBox["n", "2"]}]}], ")"}]}], "+", 
RowBox[{"dd\[Delta]", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["m", "2"], "+", 
RowBox[{
SuperscriptBox["l", "2"], " ", 
SuperscriptBox["n", "2"]}]}], ")"}]}]}], 
RowBox[{
RowBox[{
FractionBox["3", "2"], " ", "dd\[Sigma]", " ", "l", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n"}], "+", 
RowBox[{"dd\[Pi]", " ", "l", " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}]}]}], ")"}], " ", "n"}], "-", 
RowBox[{"dd\[Delta]", " ", "l", " ", 
RowBox[{"(", 
RowBox[{"1", "+", 
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "+", 
SuperscriptBox["m", "2"]}], ")"}]}]}], ")"}], " ", "n"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", 
FractionBox["1", "2"]}], " ", 
SqrtBox["3"], " ", "dd\[Delta]", " ", "l", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n"}], "+", 
RowBox[{
SqrtBox["3"], " ", "dd\[Pi]", " ", "l", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"], "-", 
SuperscriptBox["n", "2"]}], ")"}]}], "+", 
RowBox[{
SqrtBox["3"], " ", "dd\[Sigma]", " ", "l", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}]}]}]},
{
RowBox[{
FractionBox["1", "2"], " ", 
SqrtBox["3"], " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "sd\[Sigma]R"}], 
RowBox[{
RowBox[{
RowBox[{"-", "l"}], " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"]}], ")"}], " ", "pd\[Pi]R"}], "-", 
RowBox[{
FractionBox["1", "2"], " ", 
SqrtBox["3"], " ", "l", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{"m", " ", 
RowBox[{"(", 
RowBox[{"1", "+", 
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "pd\[Pi]R"}], "-", 
RowBox[{
FractionBox["1", "2"], " ", 
SqrtBox["3"], " ", "m", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n", " ", "pd\[Pi]R"}], "-", 
RowBox[{
FractionBox["1", "2"], " ", 
SqrtBox["3"], " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n", " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", "dd\[Delta]", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
RowBox[{
FractionBox["3", "2"], " ", "dd\[Sigma]R", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
RowBox[{"2", " ", "dd\[Pi]R", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "+", 
SuperscriptBox["m", "2"]}], ")"}]}]}], 
RowBox[{
RowBox[{
FractionBox["3", "2"], " ", "dd\[Sigma]R", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n"}], "+", 
RowBox[{"dd\[Delta]", " ", "m", " ", 
RowBox[{"(", 
RowBox[{"1", "+", 
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}]}]}], ")"}], " ", "n"}], "-", 
RowBox[{"dd\[Pi]R", " ", "m", " ", 
RowBox[{"(", 
RowBox[{"1", "+", 
RowBox[{"2", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}]}]}], ")"}], " ", "n"}]}], 
RowBox[{
RowBox[{
FractionBox["3", "2"], " ", "dd\[Sigma]R", " ", "l", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n"}], "+", 
RowBox[{"dd\[Pi]R", " ", "l", " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{"2", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}]}]}], ")"}], " ", "n"}], "-", 
RowBox[{"dd\[Delta]", " ", "l", " ", 
RowBox[{"(", 
RowBox[{"1", "+", 
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "+", 
SuperscriptBox["m", "2"]}], ")"}]}]}], ")"}], " ", "n"}]}], 
RowBox[{
RowBox[{
FractionBox["3", "4"], " ", "dd\[Sigma]", " ", 
SuperscriptBox[
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], "2"]}], "+", 
RowBox[{"dd\[Pi]", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"], "-", 
SuperscriptBox[
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], "2"]}], ")"}]}], "+", 
RowBox[{"dd\[Delta]", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "4"], " ", 
SuperscriptBox[
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], "2"]}], "+", 
SuperscriptBox["n", "2"]}], ")"}]}]}], 
RowBox[{
RowBox[{
SqrtBox["3"], " ", "dd\[Pi]", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "+", 
SuperscriptBox["m", "2"]}], ")"}], " ", 
SuperscriptBox["n", "2"]}], "+", 
RowBox[{
FractionBox["1", "4"], " ", 
SqrtBox["3"], " ", "dd\[Delta]", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", 
RowBox[{"(", 
RowBox[{"1", "+", 
SuperscriptBox["n", "2"]}], ")"}]}], "+", 
RowBox[{
FractionBox["1", "2"], " ", 
SqrtBox["3"], " ", "dd\[Sigma]", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}]}]}]},
{
RowBox[{
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}], " ", "sd\[Sigma]R"}], 
RowBox[{
RowBox[{
SqrtBox["3"], " ", "l", " ", 
SuperscriptBox["n", "2"], " ", "pd\[Pi]R"}], "-", 
RowBox[{"l", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}], " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{
SqrtBox["3"], " ", "m", " ", 
SuperscriptBox["n", "2"], " ", "pd\[Pi]R"}], "-", 
RowBox[{"m", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}], " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", 
SqrtBox["3"]}], " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n", " ", "pd\[Pi]R"}], "-", 
RowBox[{"n", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}], " ", "pd\[Sigma]R"}]}], 
RowBox[{
RowBox[{
RowBox[{"-", "2"}], " ", 
SqrtBox["3"], " ", "dd\[Pi]R", " ", "l", " ", "m", " ", 
SuperscriptBox["n", "2"]}], "+", 
RowBox[{
FractionBox["1", "2"], " ", 
SqrtBox["3"], " ", "dd\[Delta]", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{"1", "+", 
SuperscriptBox["n", "2"]}], ")"}]}], "+", 
RowBox[{
SqrtBox["3"], " ", "dd\[Sigma]R", " ", "l", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}]}]}], 
RowBox[{
RowBox[{
RowBox[{"-", 
FractionBox["1", "2"]}], " ", 
SqrtBox["3"], " ", "dd\[Delta]", " ", "m", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n"}], "+", 
RowBox[{
SqrtBox["3"], " ", "dd\[Pi]R", " ", "m", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"], "-", 
SuperscriptBox["n", "2"]}], ")"}]}], "+", 
RowBox[{
SqrtBox["3"], " ", "dd\[Sigma]R", " ", "m", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}]}]}], 
RowBox[{
RowBox[{
RowBox[{"-", 
FractionBox["1", "2"]}], " ", 
SqrtBox["3"], " ", "dd\[Delta]", " ", "l", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"]}], ")"}], " ", "n"}], "+", 
RowBox[{
SqrtBox["3"], " ", "dd\[Pi]R", " ", "l", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"], "-", 
SuperscriptBox["n", "2"]}], ")"}]}], "+", 
RowBox[{
SqrtBox["3"], " ", "dd\[Sigma]R", " ", "l", " ", "n", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}]}]}], 
RowBox[{
RowBox[{
SqrtBox["3"], " ", "dd\[Pi]R", " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "+", 
SuperscriptBox["m", "2"]}], ")"}], " ", 
SuperscriptBox["n", "2"]}], "+", 
RowBox[{
FractionBox["1", "4"], " ", 
SqrtBox["3"], " ", "dd\[Delta]", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", 
RowBox[{"(", 
RowBox[{"1", "+", 
SuperscriptBox["n", "2"]}], ")"}]}], "+", 
RowBox[{
FractionBox["1", "2"], " ", 
SqrtBox["3"], " ", "dd\[Sigma]R", " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "-", 
SuperscriptBox["m", "2"]}], ")"}], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
RowBox[{"-", 
SuperscriptBox["l", "2"]}], "-", 
SuperscriptBox["m", "2"]}], ")"}]}], "+", 
SuperscriptBox["n", "2"]}], ")"}]}]}], 
RowBox[{
RowBox[{
SuperscriptBox[
RowBox[{"(", 
RowBox[{
SuperscriptBox["n", "2"], "-", 
RowBox[{
FractionBox["1", "2"], " ", 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"]}], ")"}]}]}], ")"}], "2"], "dd\[Sigma]"}], " ", "+", 
RowBox[{"3", " ", 
SuperscriptBox["n", "2"], 
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"]}], ")"}], " ", "dd\[Pi]"}], " ", "+", 
RowBox[{
FractionBox["3", "4"], "  ", 
SuperscriptBox[
RowBox[{"(", 
RowBox[{
SuperscriptBox["l", "2"], "+", 
SuperscriptBox["m", "2"]}], ")"}], "2"], "dd\[Delta]"}]}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\);
];
];
];
(*all real, transpose is ok*)
(*pureHint=RU+RU\[Transpose]-RU IdentityMatrix[9]//N;*)
pureHint[[1;;MAXORB[a1],1;;MAXORB[a2]]]//N
]


(* ::Input::Initialization:: *)
SOC[i_,\[Xi]p_,\[Xi]sp_]:=Module[{a=types[[i]],s0,sx,sy,sz,ht},
s0={ {0, 0}, {0, 0}};
sx={ {0, 1},{1, 0}};
sy={{0, -I}, {I, 0}};
sz={{1, 0},{0, -1}};
ht={
{s0, s0, s0, s0},
{s0, s0, -I sz \[Xi]p, I sy \[Xi]p},
{s0, I sz \[Xi]p, s0, -I sx \[Xi]p},
{s0, - I sy \[Xi]p, I sx \[Xi]p, s0}
}//ArrayFlatten;
ht[[1;;MAXORB[a]*2,1;;MAXORB[a]*2]]
];

DFTBSolver[atoms_,types_,bonds_,spinOn:False_]:=Module[{natom,h},
natom=Length[atoms];
(*coord[i_]:=Module[{x,y,z},
{x,y,z}=atoms[[i]];
{x,y,z}=Normalize[{x,y,0}];
({
 {1, 0, 0},
 {0, 1, 0},
 {0, 0, 1}
})
];*)
R[i_,j_]:=atoms[[i]]-atoms[[j]];
(*\:5bfc\:5165\:6240\:6709hopping\:77e9\:9635\:5143*)
(*hopping H*)

Hint[i_,j_]:=Module[{a1,a2},
a1=types[[i]];a2=types[[j]];
Slater[R[j,i],a1,a2]
];
(*onsite H*)
Expand22[m_]:=Map[{{1,0},{0,1}}#&,m,{2}]//ArrayFlatten;

Hon[i_]:=Module[{a=types[[i]]},
DiagonalMatrix[
{
eon[a~~"."~~"s"],
eon[a~~"."~~"p"],eon[a~~"."~~"p"],eon[a~~"."~~"p"],
eon[a~~"."~~"d"],eon[a~~"."~~"d"],eon[a~~"."~~"d"],eon[a~~"."~~"d"],eon[a~~"."~~"d"]
},0,MAXORB[a]
]//N
];








Hempty[i_,j_]:=Array[0#&,{MAXORB[types[[i]]],MAXORB[types[[j]]]}];

HSO[i_,j_]:=\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
RowBox[{
RowBox[{"Expand22", "[", 
RowBox[{"Hon", "[", "i", "]"}], "]"}], "+", 
RowBox[{"SOC", "[", 
RowBox[{"i", ",", "\[Xi]p", ",", "\[Xi]sp"}], "]"}]}], 
RowBox[{"i", "==", "j"}]},
{
RowBox[{
RowBox[{"Hint", "[", 
RowBox[{"i", ",", "j"}], "]"}], "//", "Expand22"}], 
RowBox[{
RowBox[{"MemberQ", "[", 
RowBox[{"bonds", ",", 
RowBox[{"{", 
RowBox[{"i", ",", "j"}], "}"}]}], "]"}], "||", 
RowBox[{"MemberQ", "[", 
RowBox[{"bonds", ",", 
RowBox[{"{", 
RowBox[{"j", ",", "i"}], "}"}]}], "]"}]}]},
{
RowBox[{
RowBox[{"Hempty", "[", 
RowBox[{"i", ",", "j"}], "]"}], "//", "Expand22"}], "True"}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\);
H[i_,j_]:=\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
RowBox[{"Hon", "[", "i", "]"}], 
RowBox[{"i", "==", "j"}]},
{
RowBox[{"Hint", "[", 
RowBox[{"i", ",", "j"}], "]"}], 
RowBox[{
RowBox[{"MemberQ", "[", 
RowBox[{"bonds", ",", 
RowBox[{"{", 
RowBox[{"i", ",", "j"}], "}"}]}], "]"}], "||", 
RowBox[{"MemberQ", "[", 
RowBox[{"bonds", ",", 
RowBox[{"{", 
RowBox[{"j", ",", "i"}], "}"}]}], "]"}]}]},
{
RowBox[{"Hempty", "[", 
RowBox[{"i", ",", "j"}], "]"}], "True"}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\);
(*Execution*)
(*Print[natom];*)
(*Print[spinOn];*)
If[spinOn,
h=Array[HSO,{natom,natom}]//ArrayFlatten;
h,
h=Array[H,{natom,natom}]//ArrayFlatten
];
(*Print[Dimensions[h]];*)
SparseArray[h]
]

DFTBSolver[atoms_,types_]:=DFTBSolver[atoms,types,GetBonds[atoms,types],False];


DFTBSolver2[atoms_,types_,target_:"H"]:=Module[{
atoms1,atoms2,
bonds,
natom,
orbInAtom,
norb,
orbInForAtomIn,
HintBands,HonBands,
SintBands,SonBands,
h,s,
(*private function*)
R,Hint,Hon,Sint,Son
},
(*if two sets of atom coordinates get in, use 1st for bra and 2nd for ket*)
If[Length[Dimensions[atoms]]==3,
atoms1=atoms[[1]];atoms2=atoms[[2]],
atoms1=atoms2=atoms];

bonds=GetBonds[atoms1,types];
natom=Length[atoms1];
orbInAtom=MAXORB/@types;
norb=Total[MAXORB/@types];

R[i_,j_]:=atoms2[[i]]-atoms1[[j]];
(*\:5bfc\:5165\:6240\:6709hopping\:77e9\:9635\:5143*)
Hint[i_,j_]:=Slater[R[j,i],types[[i]],types[[j]]];
Sint[i_,j_]:=SlaterOver[R[j,i],types[[i]],types[[j]]];

Hon[i_]:=Module[{a=types[[i]]},
DiagonalMatrix[
{
eon[a~~"."~~"s"],
eon[a~~"."~~"p"],eon[a~~"."~~"p"],eon[a~~"."~~"p"],
eon[a~~"."~~"d"],eon[a~~"."~~"d"],eon[a~~"."~~"d"],eon[a~~"."~~"d"],eon[a~~"."~~"d"]
},0,MAXORB[a]
]//N
];
Son[i_]:=IdentityMatrix[MAXORB[types[[i]]]];
orbInForAtomIn=Accumulate[{1}~Join~orbInAtom][[1;;-2]];
If[target=="H",
HintBands=Band[{orbInForAtomIn[[#[[1]]]],orbInForAtomIn[[#[[2]]]]}]->Hint[#[[1]],#[[2]]]&/@bonds;
HonBands=Band[{orbInForAtomIn[[#]],orbInForAtomIn[[#]]}]->Hon[#]&/@Range[natom];
h=SparseArray[Join[HintBands,HonBands],{norb,norb}],

SintBands=Band[{orbInForAtomIn[[#[[1]]]],orbInForAtomIn[[#[[2]]]]}]->Sint[#[[1]],#[[2]]]&/@bonds;
SonBands=Band[{orbInForAtomIn[[#]],orbInForAtomIn[[#]]}]->Son[#]&/@Range[natom];
s=SparseArray[Join[SintBands,SonBands],{norb,norb}]
]

]



DFTBSolverOver[atoms_,types_]:=Module[{bonds,natom,h},
natom=Length[atoms];
bonds=GetBonds[atoms,types];

R[i_,j_]:=atoms[[i]]-atoms[[j]];

Son[i_]:=IdentityMatrix[MAXORB[types[[i]]]];
Sint[i_,j_]:=Module[{a1,a2},
a1=types[[i]];a2=types[[j]];
SlaterOver[R[j,i],a1,a2]
];
Sempty[i_,j_]:=Array[0#&,{MAXORB[types[[i]]],MAXORB[types[[j]]]}];

OL[i_,j_]:=\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
RowBox[{"Son", "[", "i", "]"}], 
RowBox[{"i", "==", "j"}]},
{
RowBox[{"Sint", "[", 
RowBox[{"i", ",", "j"}], "]"}], 
RowBox[{
RowBox[{"MemberQ", "[", 
RowBox[{"bonds", ",", 
RowBox[{"{", 
RowBox[{"i", ",", "j"}], "}"}]}], "]"}], "||", 
RowBox[{"MemberQ", "[", 
RowBox[{"bonds", ",", 
RowBox[{"{", 
RowBox[{"j", ",", "i"}], "}"}]}], "]"}]}]},
{
RowBox[{"Sempty", "[", 
RowBox[{"i", ",", "j"}], "]"}], "True"}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\);

s=Array[OL,{natom,natom}]//ArrayFlatten;
N[SparseArray[s]]
]


(* ::Section:: *)
(*Plot DFTB result*)


PlotEigen[eval_]:=Module[{listP},
listP=Table[{i,eval[[i]]},{i,Length[eval]}];
ListPlot[listP]
]

(*wave function plot*)
wfp[array_,pos_,\[Mu]_]:=Module[{arrayT},
arrayT=array~ArrayReshape~{NATOM,4};
Graphics3D[
{
Table[{
Text[Style[i,Black,Thick,16],pos[[i]]],
Sphere[pos[[i]],arrayT[[i,\[Mu]]]*4]
},
{i,Length[pos]}
]
},
ViewPoint->Front,
ImageSize->200
]
]

CalcGap[eval_,NValence_]:=Module[{},eval[[NValence+1]]-eval[[NValence]]]
CalcHLAV[eval_,NValence_]:=Module[{},(eval[[NValence]]+eval[[NValence+1]])/2]



CalcV12[{atoms1_,types1_},{atoms2_,types2_}]:=Module[{s1,s2,atoms,types,bonds,h},
s1=Total[MAXORB/@types1];
s2=Total[MAXORB/@types2];
atoms=Join[atoms1,atoms2];
types=Join[types1,types2];
bonds=GetBonds[atoms,types];
h=DFTBSolver2[atoms,types];
h[[1;;s1,-s2;;-1]]
]

CalcS12[{atoms1_,types1_},{atoms2_,types2_}]:=Module[{s1,s2,atoms,types,bonds,s},
s1=Total[MAXORB/@types1];
s2=Total[MAXORB/@types2];
atoms=Join[atoms1,atoms2];
types=Join[types1,types2];
bonds=GetBonds[atoms,types];
s=DFTBSolverOver[atoms,types];
s[[1;;s1,-s2;;-1]]
]


(* ::Section:: *)
(*Repulsive Potential*)


(*\:5bfc\:5165\:6240\:6709\:6392\:65a5\:52bf*)
GetREP[bond_,prefix_]:=Module[{atom1,atom2,NSP,text,nInt,spline,it},
path=prefix~~bond~~".skf";
text=Import[path,"Lines"];
(*Assert file name must be A-B*)
{atom1,atom2}=bond~StringSplit~"-";
NSP=Position[text,"Spline"]//First//First;
{nInt,cutoff}=text[[NSP+1]]//StringSplit//ToExpression;
{a1,a2,a3}=text[[NSP+2]]//StringSplit//ToExpression;
spline=text[[NSP+3;;NSP+nInt+2]]//StringSplit//ToExpression;
lsp1=spline[[1,1]];
it=Interpolation[spline[[All,{1,3}]],InterpolationOrder->1];
(*\:8ddd\:79bb\:5c0f\:4e8e\:7b2c\:4e00\:53d6\:6837\:957f\:5ea6\:7528\:6307\:6570\:8868\:8fbe\:5f0f\:ff0c\:4ecb\:4e8e\:7b2c\:4e00\:53d6\:6837\:957f\:5ea6\:548c\:622a\:65ad\:957f\:5ea6\:4e4b\:95f4\:7528\:7ebf\:6027\:63d2\:503c*)
REP[bond,r_]=Module[{},s=r/AUL;f=\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
RowBox[{"Exp", "[", 
RowBox[{
RowBox[{
RowBox[{"-", "a1"}], "*", "s"}], "+", "a2"}], "]"}], 
RowBox[{"s", "<", "lsp1"}]},
{
RowBox[{"it", "[", "s", "]"}], 
RowBox[{"lsp1", "<=", "s", "<", 
RowBox[{"cutoff", "-", "lsp1"}]}]},
{"0", "True"}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)];f*AUE(*\:91cd\:6574\:5355\:4f4d*)
];
(*\:8ba1\:7b97\:539f\:5b50\:4e4b\:95f4\:6392\:65a5\:52bf*)
(*HREP[i_,j_]:=Module[{r,bond,rep},
r=R[i,j]//Norm;
bond=types[[i]]~~"-"~~types[[j]];
rep=REP[bond,r];
Table[rep,{4},{4}]
];*)


GetRepArr[atoms_,types_]:=Block[{x,y,z,dx,dy,dz,dr,rep},
{x,y,z}=atoms\[Transpose];
{dx,dy,dz}=Outer[Subtract,#,#]&/@{x,y,z};
dr=Sqrt[dx^2+dy^2+dz^2];
rep=Table[
If[0<dr[[i,j]]<4,
REP[types[[i]]~~"-"~~types[[j]],dr[[i,j]]],
0
],
{i,Length[dr]},{j,Length[dr]}
]//SparseArray;
rep
]

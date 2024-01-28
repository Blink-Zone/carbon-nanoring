(* ::Package:: *)

(* ::Section:: *)
(*XYZT*)


TranslateXYZT[xyzt_,v_]:=Append[xyzt[[1;;3]]+v,xyzt[[4]]]


ImportXYZFrame[fileXYZ_]:=Block[{string},
(*\:5bfc\:5165xyz\:5e27\:5e8f\:5217*)
fullTable=Import[fileXYZ,"Table"];
atomNum=fullTable[[1,1]];
frames=fullTable~Partition~(atomNum+2);
atomsFrame=frames[[All,3;;-1,2;;4]];
typesFrame=frames[[All,3;;-1,1]];
{atomsFrame,typesFrame}
]
ExportXYZFrame[path_,types_,atomsFrame_]:=Block[{string},
string=({ExportString[{types,#},{"XYZ",{"VertexTypes","VertexCoordinates"}}],"\r"}&/@atomsFrame)//StringJoin;
Export[path,string,"Text"]
]


SortXYZFrameByCertainZ[in_,out_,n_]:=Module[{},
frames=ImportXYZFrame[in];
framesT=Transpose[frames,{3,2,1}];
framesTS=SortBy[framesT,#[[n,1,3]]&];
framesTST=Transpose[framesTS,{3,2,1}];
ExportXYZFrame[out,framesTST[[2,1]],framesTST[[1]]];
]


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
table=SortBy[table,#[[2;;4]]&];
types=table[[All,1]];
atoms=table[[All,2;;4]];
CopyFile[f,f~~".bk"];
ExportXYZ[f,atoms,types]
]

MOL2XYZ[fileMOL_,fileXYZ_]:=Module[{atoms,types,bonds},
{atoms,types,bonds}=ImportMOL[fileMOL];
ExportXYZ[fileXYZ,atoms,types]
]


RotateSystem[atomsA_,\[Theta]_,range_]:=Module[{atomsB,atomsC,axis={0,0,0}},
atomsC=atomsA;
atomsB=RotationTransform[\[Theta],{0,0,1},axis]/@atomsA;
atomsC[[range[[1]];;range[[2]]]]=atomsB[[range[[1]];;range[[2]]]];
atomsC
]

RotateSystem2[atomsA_,rotRange_,rotAxis_,rotAngle_,rotPivot_:{0,0,0}]:=Module[{atomsB,atomsC},
atomsC=atomsA;
atomsB=RotationTransform[rotAngle,rotAxis,rotPivot]/@atomsA;
atomsC[[rotRange[[1]];;rotRange[[2]]]]=atomsB[[rotRange[[1]];;rotRange[[2]]]];
atomsC
]


hrot[\[Theta]_]:=hrot[\[Theta]]=DFTBSolver[RotateSystem2[atoms,{21,-21},r2-r1,\[Theta],r1],types];
mrot[\[Theta]_]:=mrot[\[Theta]]=(hrot[\[Theta]+d\[Theta]]-hrot[\[Theta]-d\[Theta]])/(2d\[Theta])

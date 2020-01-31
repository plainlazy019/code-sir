function [] = FEM_Class_1D();

% uses \ to get inverse 
% arranges as NQ * NQ matrix

% ONE degree of freedom per node.
close all
clear all
clc

InputData;
Stiffness;
ModifyForBCon;
StressCalc;
ReactionCal;
Output;
%------------------------  function InputData  ---------------------------
function [] = InputData();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT DT S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
global NQ EL

FILE1 = input('Input Data File Name ','s');
LINP  = fopen(FILE1,'r');
FILE2 = input('Output Data File Name ','s');
LOUT  = fopen(FILE2,'w');

TITLE = fgets(LINP)
DUMMY = fgets(LINP)
TMP = str2num(fgets(LINP))
[NN, NE, NM, NDIM, NEN, NDN] = deal(TMP(1),TMP(2),TMP(3),TMP(4),TMP(5),TMP(6));

NQ = NDN * NN;

DUMMY = fgets(LINP);
TMP = str2num(fgets(LINP));
[ND, NL, NPR, EL]= deal(TMP(1),TMP(2), TMP(3), TMP(4));

%----- Coordinates -----
DUMMY = fgets(LINP);
for I=1:NN
   TMP = str2num(fgets(LINP));
   [N, X(N,:)]=deal(TMP(1),TMP(2:1+NDIM));
end
%----- Connectivity -----
DUMMY = fgets(LINP);
for I=1:NE
   TMP = str2num(fgets(LINP));
   [N,NOC(N,:), MAT(N,:), AREA(N,:)] = ...
   deal(TMP(1),TMP(2:1+NEN), TMP(2+NEN), TMP(3+NEN));
end

%----- Specified Displacements -----
DUMMY = fgets(LINP);
for I=1:ND
   TMP = str2num(fgets(LINP));
   [NU(I,:),U(I,:)] = deal(TMP(1), TMP(2));
end
%----- Component Loads -----
DUMMY = fgets(LINP);
F = zeros(NQ,1);
for I=1:NL
   TMP = str2num(fgets(LINP));
   [N,F(N)]=deal(TMP(1),TMP(2));
end

%----- Material Properties -----
DUMMY = fgets(LINP);
for I=1:NM
   TMP = str2num(fgets(LINP));
   [N, PM(N,:)] = deal(TMP(1), TMP(2:NPR+1));
end

%------------------------  function Stiffness  ---------------------------
function []=Stiffness();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT DT S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
global NQ EL

%----- First initialize stiffness matrix
S = zeros(NQ,NQ);      

%-----  Stiffness Matrix modified -----
for N = 1:NE
   p = NOC(N, 1);
   s = NOC(N, 2);
   I3 = MAT(N);

  
   EAL = PM(I3, 1) * AREA(N) / EL;

%----------- Element Stiffness Matrix SE() -----------

   SE=[1 -1; -1 1];
  

   disp(sprintf('..... Adding %dth Element Stiffness to Global Locations',N));


   row=0;
   for I=[p,s]
       col=0;
       row=row+1;
       for J=[p,s];
           col=col+1;
           S(I,J)=S(I,J)+ (SE(row,col)*EAL);
       end
   end
end
Stiffness_matrix=S


%------------------------  function ModifyForBCon  ---------------------------
function []=ModifyForBCon();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT DT S
global PM NU U MPC BT STRESS REACT
global CNST J
global TITLE FILE1 FILE2 FF
global LINP LOUT
global NQ



%----- Modify for Boundary Conditions -----
%    --- Displacement BC ---
BCC=[];
SS=S;
FF=F;
FO=F
v=[];
BCC=zeros(NQ,1)
for I = 1:ND
    NI = NU(I);
    v=[v,NI];
    BC=S(:,NI);
    BCC=BCC+(BC*U(I));
end
 
u=[];
for I=1:NQ
    if I~=v
        u=[u I];
    end
end

FF=F(u);
FF=FF-BCC(u);
F(v)=U;
SS(v,:)=[];
SS(:,v)=[];
[FF]=SS\FF
F(u)=FF

%------------------------  function StressCalc  ---------------------------
function []=StressCalc();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT DT S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
global NQ EL


%----- Stress Calculation modified for one degree of freedom system -----

for I = 1:NE
   I1 = NOC(I, 1);
   I2 = NOC(I, 2);
   I3 = MAT(I);
 
   J1 = I1;
   K1 = I2;
   DLT = F(K1) - F(J1);
   STRESS(I) = PM(I3, 1) * (DLT / EL);
end 



%------------------------  function ReactionCal  ---------------------------
function []=ReactionCal();
global NN NE NM NDIM NEN NDN
global ND NL NCH NPR NMPC NBW
global X NOC F AREA MAT DT S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT
global NQ

%----- Reaction Calculation -----
for I = 1:ND
    NI = NU(I);    
    BC=S(:,NI)
   REACT(I)=BC'*F;
end

%------------------------  function Output  ---------------------------
function []=Output();
global NN NE NM NDIM NEN NDN
global ND NL NMPC NBW
global X NOC F AREA MAT DT S
global PM NU U MPC BT STRESS REACT
global CNST
global TITLE FILE1 FILE2
global LINP LOUT

disp(TITLE);
fprintf(LOUT,'%s\n',TITLE);
disp('NODE#   DISPLACEMENT');
fprintf(LOUT,'NODE#   DISPLACEMENT\n');
for I=1:NN
	disp(sprintf('%d %15.5G', I, F(I)));
	fprintf(LOUT,'%d %15.5G\n', I, F(I));
end
%----- Stress  -----
disp('ELEM#    STRESS');
fprintf(LOUT,'ELEM#    STRESS\n');
for N=1:NE
	disp(sprintf('%d %15.5G',N, STRESS(N)));
	fprintf(LOUT,'%d %15.5G\n',N, STRESS(N));
end
%----- Reaction  -----
disp('NODE#    REACTION');
fprintf(LOUT,'NODE#    REACTION\n');
for I=1:ND
   N = NU(I);
	disp(sprintf('%d %15.5G',N, REACT(I)));
	fprintf(LOUT,'%d %15.5G\n', N, REACT(I));
end

fclose(LOUT);
disp(sprintf('RESULTS ARE IN FILE : %s', FILE2));
% Truss2D_Ex02
%----------------------------------------------------------------
% PURPOSE 
%    Static Non-Linear analysis o a 2D steel truss.
%
%----------------------------------------------------------------
%
% LAST MODIFIED: L.F.Veduzco    2023-05-25
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------

clc
clear all

%% Materials
E=[29000*1000*0.4536/2.54^2; % Elements' Modulus of Elasticity
    29000*1000*0.4536/2.54^2;
    29000*1000*0.4536/2.54^2;
   29000*1000*0.4536/2.54^2;
   29000*1000*0.4536/2.54^2;
   29000*1000*0.4536/2.54^2;
   29000*1000*0.4536/2.54^2;
   29000*1000*0.4536/2.54^2;
   29000*1000*0.4536/2.54^2;
   29000*1000*0.4536/2.54^2;
   29000*1000*0.4536/2.54^2];

Fy=[3515; % Yield stress
    3515;
    3515;
    3515;
    3515;
    3515;
    3515;
    3515;
    3515;
    3515;
    3515];

%% Geometry
area1=30.26; r1=5.16; % OC 152 x 7.1
area2=50.71; r2=5.59; % OC 168 x 11.0

A=[area2; % element's cross-section area
    area1;
    area1;
    area2;
    area2;
    area2;
    area2;
    area2;
    area2;
    area1;
    area1];

r=[r2; % element's cross-section radius rotation
   r1;
   r1;
   r2;
   r2;
   r2;
   r2;
   r2;
   r2;
   r1;
   r1];

%% Node coordinates
coordxy=[0,0;
         300,0;
         600,0;
         900,0;
         150,400;
         450,800;
         750,400]; 

%% Topology
% Initial and final node of each bar
ni=[1;2;3;1;5;2;6;3;7;5;6];
nf=[2;3;4;5;2;6;3;7;4;6;7];

%% Boundary conditions
bc=[1 0;
    2 0;
    7 0;
    8 0];

%% Loads
load=200; %kgf
load01=2*load;
load02=-2*load;
load03=load;
load04=-2*load;

initialLoads=[load01; load02; load04];
edofLoads=[9;10;12];

%% Pushover analysis
[finalLoads,collecDxNodes,collecDyNodes,collecFS,Uglobal,criticalBarsSteps,...
 critical_ten_comp_bars_global,dam,defBarHistory,StrBarHistory]=...
 ElastoPlasticPushover2DTruss(E,Fy,A,r,coordxy,ni,nf,bc,edofLoads,...
 initialLoads,10,1,[1,9,10],0.001)

%% Unloading - Residual stresses and strains
[residStr,residDef,defBarHist,StrBarHist]=Unload2DTrusses...
(defBarHistory,StrBarHistory,finalLoads,edofLoads,[1,9,10],E,A,...
coordxy,ni,nf,bc)

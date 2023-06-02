% Truss2D_Ex01
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

clear all
clc

%% Materials
E=[29000*1000*0.4536/2.54^2; % Elements' Modulus of Elasticity (Kg/cm2)
    29000*1000*0.4536/2.54^2;
    29000*1000*0.4536/2.54^2;
   29000*1000*0.4536/2.54^2;
   29000*1000*0.4536/2.54^2;
   29000*1000*0.4536/2.54^2;
   29000*1000*0.4536/2.54^2;
   29000*1000*0.4536/2.54^2;
   29000*1000*0.4536/2.54^2];

Fy=[3515; % Yield stress (Kg/cm2)
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

A=[area2; % element's cross-section area (cm2)
    area1;
    area1;
    area1;
    area1;
    area2;
    area1;
    area1;
    area2];

r=[r2; % element's cross-section radius rotation (cm)
   r1;
   r1;
   r1;
   r1;
   r2;
   r1;
   r1;
   r2];

%% Topology
% Node coordinates (cm)
coordxy=[0,0;
         600,0;
         1100,0;
         0,600;
         600,500;
         1100,500]; 

% Initial and final node of each bar
ni=[1;4;4;5;5;2;2;6;1];
nf=[2;2;5;2;6;6;3;3;5];

%% Boundary conditions
bc=[1 0;
    2 0;
    7 0;
    8 0];

%% Loads
load=200; %kgf 
load01=load;
load02=2*load;
load03=2*load;
        
initialLoads=[load01; load02; load03];
edofLoads=[4;6;11];

%% Pushover analysis
[finalLoads,collecDxNodes,collecDyNodes,collecFS,Uglobal,criticalBarsSteps,...
 critical_ten_comp_bars_global,dam,defBarHistory,StrBarHistory]=...
 ElastoPlasticPushover2DTruss(E,Fy,A,r,coordxy,ni,nf,...
 bc,edofLoads,initialLoads,10,1,[1,2,3],0.001)

%% Unloading - Residual stresses and strains
[residStr,residDef,defBarHist,StrBarHist]=Unload2DTrusses...
(defBarHistory,StrBarHistory,finalLoads,edofLoads,[1,2,3],E,A,...
coordxy,ni,nf,bc)
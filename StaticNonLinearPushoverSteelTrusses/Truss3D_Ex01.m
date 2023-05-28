% Truss3D_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    Static Non-Linear analysis o a 3D steel truss.
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
E=[29000*1000*0.4536/2.54^2; % Elements' Modulus of Elasticity
    29000*1000*0.4536/2.54^2;
    29000*1000*0.4536/2.54^2;
   29000*1000*0.4536/2.54^2];

Fy=[3515; % Yield stress
    3515;
    3515;
    3515];

%% Geometry
area1=30.26; r1=5.16; % OC 152 x 7.1 
area2=50.71; r2=5.59; % OC 168 x 11.0

A=[area2; % element's cross-section area
    area1;
    area2;
    area1];

r=[r2; % element's cross-section radius rotation
   r1;
   r1;
   r2];

%% Topology
% Node coordinates
coordxyz=[0,0,500;
         -400,-200,0;
         0,400,0;
         400,400,0;
         400,-400,0]; 

% Initial and final node of each bar
ni=[1;1;1;1];
nf=[2;3;4;5];

%% Boundary conditions
bc=[4 0;
    5 0;
    6 0;
    7 0;
    8 0;
    9 0;
    10 0;
    11 0;
    12 0;
    13 0;
    14 0;
    15 0];

%% Loads
load=200; %kgf
load01=load;
load02=-2*load;
        
initialLoads=[load01; load02];
edofLoads=[1;3];

%% Pushover analysis
[collecDxNodes,collecDyNodes,collecDzNodes,collectionFS,Uglobal,...
 criticalBarsGlobal,criticalTenCompBarsGlobal,dam,defBarHistory,...
 strBarHistory]=ElastoPlasticPushover3DTruss(E,Fy,A,r,coordxyz,ni,nf,...
 bc,edofLoads,initialLoads,10,1,[1,2,3,4],0.001)

function [residStr,residDef,defBarHistUnl,StrBarHistUnl]=Unload3DTrusses...
    (defBarHist,StrBarHist,finalLoads,edofLoads,barsplotPD,E,A,...
    coordxyz,ni,nf,bc)
%------------------------------------------------------------------------
% Syntax:
% [residStr,residDef,defBarHist,StrBarHist]=Unload3DTrusses...
%  (defBarHist,StrBarHist,finalLoads,edofLoads,barsplotPD,E,A,...
%  coordxyz,ni,nf,bc)
%
%------------------------------------------------------------------------
% PURPOSE
%  To compute the residual stress-deformations of N number of bars
%  composing a 3D truss through an unloading process from the results 
%  obtained in a previous Pushover analysis. 
%  
% 
% INPUT:  A:                     element's cross-section area
%
%         E:                     Element's Modulus of Elasticity
%
%         coordxy:               node coordinates [x,y,z]
%
%         ni,nf:                 initial and final nodes of each bar
%
%         bc:                    boundary condition vector. 
%                                Format: [n-dof, condition]
%
%         edofLoads:             list of degrees of freedom subject to
%                                the applied punctual forces
%
%         finalLoads:            punctual initial forces at each listed dof
%                                given in "edofLoads" in the same order
%
%         barsplotPD:            list of chosen bars to analyse their load
%                                history along the Pushover assessment
%
%         defBarHist:           Axial deformation history of each bar 
%                               during the pushover analysis
%
%                          [def-step-1, def-step-2, ... ]
%
%         StrBarHist:           Axial stress history of each bar 
%                               during the pushover analysis 
%
%                          [str-step-1, str-step-2, ... ]
%
% OUTPUT: defBarHistUnl:        Axial deformation history of each bar 
%                               after unloading
%
%         StrBarHistUnl:        Axial stress history of each bar 
%                               after unloading
%
%         residStr,residDef:    Residual stresses and deformation of each
%                               bar after unloading
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-06-02
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

nnodes=length(coordxyz(:,1));
nbars=length(E);

%%%--------------------------- Topology -------------------------%%%

Edof=zeros(nbars,7);
for i=1:nbars
    Edof(i,1)=i;
    Edof(i,2)=ni(i)*3-2;
    Edof(i,3)=ni(i)*3-1;
    Edof(i,4)=ni(i)*3;
    
    Edof(i,5)=nf(i)*3-2;
    Edof(i,6)=nf(i)*3-1;
    Edof(i,7)=nf(i)*3;
end
L=sqrt((coordxyz(ni,1)-coordxyz(nf,1)).^2+...
      (coordxyz(ni,2)-coordxyz(nf,2)).^2+...
      (coordxyz(ni,3)-coordxyz(nf,3)).^2); % Elements' length
  
ex=coordxyz(:,1);
ey=coordxyz(:,2);
ez=coordxyz(:,3);

f=zeros(3*nnodes,1);
f(edofLoads)=-finalLoads; % forces acting in the opposite direction

%% Static Linear Analysis
% The final loads obtained from a Pushover analysis to reach a collapse
% mechanism are used. They are applied with the same magnitude but in
% opposite direction so that a total unload is carried out

% Assembling global stiffness matrix
K=zeros(3*nnodes,3*nnodes);
for i=1:nbars
    [Ke]=bar3e([ex(ni(i)),ex(nf(i))],[ey(ni(i)),ey(nf(i))],[ez(ni(i)),...
        ez(nf(i))],[E(i),A(i)]);
    [K]=assem(Edof(i,:),K,Ke);
end
% Solution to the system of equations
[Uglobal,reaction]=solveq(K,f,bc);    

%%%%%%%-------------------Elemental forces---------------%%%%%%%%
for i=1:nbars
    ue=[Uglobal(3*ni(i)-2);
           Uglobal(3*ni(i)-1);
           Uglobal(3*ni(i));
           Uglobal(3*nf(i)-2);
           Uglobal(3*nf(i)-1);
           Uglobal(3*nf(i))];

    % Normal mechanical elements
    [esBarsNormal(:,i)]=bar3s([ex(ni(i)),ex(nf(i))],[ey(ni(i)),ey(nf(i))],...
        [ez(ni(i)),ez(nf(i))],[E(i),A(i)],ue');
    
end

%% Stresses and deformations of bars
unlStrBar=esBarsNormal(1,:)'./A; % stresses in bars
unlDefBar=esBarsNormal(1,:)'./A./E; % unitary axial deformation in bars

%% Plotting unload stresses and strains history
if isempty(barsplotPD)==0
    %%------------------- Residual stresses -----------------%%
    plasticNsteps=length(StrBarHist(1,:))-1;
    
    residStr=StrBarHist(:,plasticNsteps+1)+unlStrBar;
    residDef=defBarHist(:,plasticNsteps+1)+unlDefBar;
    
    defBarHistUnl=[defBarHist,residDef];
    StrBarHistUnl=[StrBarHist,residStr];

    ngt9=0;
    nlet9=0;
    figure(4)
    if barsplotPD(1)>9 % This condition is just to make sure that the 
                       % bars' labels are correct in the plot
        ngt9=ngt9+1;
        bariText2(ngt9,:)=strcat('Bar ',num2str(barsplotPD(1)));
        if defBarHistUnl(barsplotPD(1),plasticNsteps)==...
            StrBarHistUnl(barsplotPD(1),plasticNsteps+1) % to identify if the bar
                                                         % has yielded

            plot(defBarHistUnl(barsplotPD(1),:),StrBarHistUnl(barsplotPD(1),:),...
                'k -','LineWidth',1.8)
        else
            plot(defBarHistUnl(barsplotPD(1),[1,plasticNsteps+2]),...
                 StrBarHistUnl(barsplotPD(1),[1,plasticNsteps+2]),...
                'k -','LineWidth',1.8)
        end
        legend(bariText2(1,:))
        hold on
    else
        nlet9=nlet9+1;
        bariText(nlet9,:)=strcat('Bar ',num2str(barsplotPD(1)));
        if StrBarHistUnl(barsplotPD(1),plasticNsteps)==...
            StrBarHistUnl(barsplotPD(1),plasticNsteps+1) % to identify if the bar
                                                         % has yielded

            plot(defBarHistUnl(barsplotPD(1),:),StrBarHistUnl(barsplotPD(1),:),...
                'k -','LineWidth',1.8)
        else
            plot(defBarHistUnl(barsplotPD(1),[1,plasticNsteps+2]),...
                 StrBarHistUnl(barsplotPD(1),[1,plasticNsteps+2]),...
                'k -','LineWidth',1.8)
        end
        legend(bariText(1,:))
        hold on
    end

    for i=2:length(barsplotPD)
        if barsplotPD(i)>9 % This condition is just to make sure that the 
                           % bars' labels are correct in the plot
            ngt9=ngt9+1;
            bariText2(ngt9,:)=strcat('Bar ',num2str(barsplotPD(i)));
            if StrBarHistUnl(barsplotPD(i),plasticNsteps)==...
                StrBarHistUnl(barsplotPD(i),plasticNsteps+1) % to identify if the
                                                             % bar has yielded
                plot(defBarHistUnl(barsplotPD(i),:),...
                    StrBarHistUnl(barsplotPD(i),:),...
                    'LineWidth',1.8,'DisplayName',bariText2(ngt9,:))
            else
                plot(defBarHistUnl(barsplotPD(i),[1,plasticNsteps+2]),...
                     StrBarHistUnl(barsplotPD(i),[1,plasticNsteps+2]),...
                     'LineWidth',1.8,'DisplayName',bariText2(ngt9,:))
            end
        else
            nlet9=nlet9+1;
            bariText(nlet9,:)=strcat('Bar ',num2str(barsplotPD(i)));
            if StrBarHistUnl(barsplotPD(i),plasticNsteps)==...
                StrBarHistUnl(barsplotPD(i),plasticNsteps+1) % to identify if the
                                                             % bar has yielded
                plot(defBarHistUnl(barsplotPD(i),:),...
                    StrBarHistUnl(barsplotPD(i),:),...
                    'LineWidth',1.8,'DisplayName',bariText(nlet9,:))
            else
                plot(defBarHistUnl(barsplotPD(i),[1,plasticNsteps+2]),...
                     StrBarHistUnl(barsplotPD(i),[1,plasticNsteps+2]),...
                     'LineWidth',1.8,'DisplayName',bariText(nlet9,:))
            end
        end
    end
    xlabel('Residual deformation')
    ylabel('Residual stress')
    title('Residual stresses-deformation (Unload)')
    hold on
else
    disp('Error: no bars to assess were given !')
    residStr=[]; residDef=[]; defBarHistUnl=[]; StrBarHistUnl=[];
    return
end

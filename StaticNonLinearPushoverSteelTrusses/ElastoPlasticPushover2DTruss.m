function [collec_dx_nodes,collec_dy_nodes,collection_FS,Uglobal,...
         criticalBarsGlobal,criticalTenCompBarsGlobal,kdam,defBarHistory,...
         StrBarHistory]=ElastoPlasticPushover2DTruss(E,Fy,A,r,coordxy,ni,...
         nf,bc,edofLoads,initialLoads,dle,plotDef,barsplotPD,kdeg)
%------------------------------------------------------------------------
% Syntax:
% [collec_dx_nodes,collec_dy_nodes,collection_FS,Uglobal,...
%  criticalBarsGlobal,criticalTenCompBarsGlobal,dam,defBarHistory,...
%  StrBarHistory]=ElastoPlasticPushover2DTruss(E,Fy,A,r,coordxy,ni,nf,...
%  bc,edofLoads,initialLoads,dle,plotDef,barsplotPD,kdeg)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE
%  To perform the Static Non-Linear analysis for 2D trusses.
%  
% INPUT:  A:                     element's cross-section area
%
%         Fy:                    Element's Yield stress
%
%         E:                     Element's Modulus of Elasticity
%
%         r:                     radius of rotation
%
%         coordxy:               node coordinates [x,y]
%
%         ni,nf:                 initial and final nodes of each bar
%
%         bc:                    boundary condition vector. 
%                                Format: [n-dof, condition]
%
%         edofLoads:             list of degrees of freedom subject to
%                                the applied punctual forces
%
%         initialLoads:          punctual initial forces at each listed dof
%                                given in "edofLoads" in the same order
%
%         dle:                   load step increment
%
%         plotDef:               parameter that indicates whether or not
%                                the plot of the deformed structure before 
%                                collapse will be shown. 0 - Not to plot
%                                                        1 - To plot
%
%         barsplotPD:            list of chosen bars to analyse their load
%                                history along the Pushover assessment.
%                                When no load history analysis is required
%                                set the parameter to an empty vector [] or
%                                to 0
%
%         kdeg:                  stifness degradation of the structure in
%                                the collapse state in relation to the 
%                                initial stifness 
%
% OUTPUT: collec_dx_nodes,
%         collec_dy_nodes:      arrays containing the displacement
%                               history of each node in the x and y
%                               direction, respectively
%
%         collection_FS:        list of load increments corresponding to
%                               to the step in which each plastification 
%                               occurs
%
%         Uglobal:              is the displacement vector of each dof
%                               right before reaching the collapse 
%                               mechanism
%
%         criticalBarsGlobal:   array containing the bar ID that is
%                               plastiied at each plastification mechanism
%
%         criticalTenCompBarsGlobal:   array containing the stress state ID
%                                      at which each bar is plastified when
%                                      a plastification occurs, either:
%
%                               1 - Tension
%                               2 - Compression
% 
%         kdam:                 stiffness degreation damage
%
%         defBarHistory:        Axial deformation history of each bar 
%                               during the pushover analysis
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-05-24
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

nnodes=length(coordxy(:,1));
nbars=length(E);

%%%--------------------------- Topology -------------------------%%%

Edof=zeros(nbars,5);
for i=1:nbars
    Edof(i,1)=i;
    Edof(i,2)=ni(i)*2-1;
    Edof(i,3)=ni(i)*2;
    
    Edof(i,4)=nf(i)*2-1;
    Edof(i,5)=nf(i)*2;
end
L=sqrt((coordxy(ni,1)-coordxy(nf,1)).^2+...
      (coordxy(ni,2)-coordxy(nf,2)).^2);

[ndof,edof]=nonRestrcDof(nnodes,bc,2);

%-------------------------------------------------------------------
%               NOTA: Tension (+), Compression (-)
%-------------------------------------------------------------------

% CRITICAL TENSION LOADS
for i=1:nbars
    critical_axial_tension_bar(i,1)=Fy(i)*A(i)/1.67;
end

% CRITICAL COMPRESSION LOADS
k=1;
klr_bar=k*L./r;
for i=1:nbars
    Fe_bar(i,1)=pi^2*E(i)/(klr_bar(i)^2); % Euler stress
    Fcr_bar(i,1)=0.877*Fe_bar(i); % Critical stress
    Pn_bar(i,1)=Fcr_bar(i)*A(i); % Effective axial resistance
    critical_axial_compresion_bar(i,1)=-Pn_bar(i)/1.67; % Design resistance
end
ex=coordxy(:,1);
ey=coordxy(:,2);

initialLoads1=initialLoads; % storing initial loads 
dl=sign(initialLoads)*dle; % load increment 

collec_dx_nodes=[zeros(1,nnodes)];
collec_dy_nodes=[zeros(1,nnodes)];
collection_FS=[0];

criticalBarsGlobal=[];
criticalTenCompBarsGlobal=[];

critical_compresion_force=0; % to initialize the while loop
critical_tension_force=0;
iter=0;
f=zeros(2*nnodes,1);
collect_Kred=[];
defBarHistory=[zeros(nbars,1)]; % To store the bars' deformation history
StrBarHistory=[zeros(nbars,1)];
while (critical_tension_force==0 && critical_compresion_force==0)
    iter=iter+1;
    if iter==1 % Initial loads
        f(edofLoads)=initialLoads;
    end
    %%%%%%%%--------------------Assembly------------------%%%%%%%%%%%
    K=zeros(2*nnodes,2*nnodes);
    Kes=[];
    for i=1:nbars
        [Ke]=bar2e([ex(ni(i)),ex(nf(i))],[ey(ni(i)),ey(nf(i))],[E(i),A(i)]);
        Kes=[Kes;Ke];
        [K,f]=assem(Edof(i,:),K,Ke,f,[0;0;0;0]);
    end

    if iter==1
        %Global reduced K matrix
        globalKreduced=K(edof,edof);
        
        %Determinant of the initial sttructure's stiffness
        det_Kred=det(globalKreduced);
        collect_Kred=[collect_Kred,det_Kred];
    end
    %Solution to the system of equations
    [Uglobal,reaction]=solveq(K,f,bc);
                  
    %%%%%%%-------------------Elemental forces---------------%%%%%%%%
    for i=1:nbars
        ue=[Uglobal(2*ni(i)-1);
           Uglobal(2*ni(i));
           Uglobal(2*nf(i)-1);
           Uglobal(2*nf(i))];
       
        ke=Kes(4*i-3:4*i,:);
        
        fe=-ke*ue; % reactions in the nodes of each bar
        
        fbe_bars(:,i)=fe;
        
        % Normal mechanical elements
        [es]=bar2s([ex(ni(i)),ex(nf(i))],[ey(ni(i)),ey(nf(i))],...
            [E(i),A(i)],ue');
        es_bars_normal(:,i)=es;
    end
    
    % STOP WHEN A MEMBER RESISTANCE IN TENSION IS OVERPAST
    critical_bars=[];
    critical_bars_step=zeros(1,nbars);
    
    critical_ten_comp_bars=zeros(1,nbars);
    critical_ten_comp=[];
    for j=1:nbars
        axialforce=es_bars_normal(1,j);
        if axialforce>0 % for axial forces in tension
            if(axialforce>critical_axial_tension_bar(j))
                critical_bars=[critical_bars,j];
                critical_bars_step(j)=j;
                critical_ten_comp=[critical_ten_comp,1];
                critical_ten_comp_bars(j)=1;
            end
        else           % for axial forces in compression
            if(axialforce<critical_axial_compresion_bar(j))
                critical_bars=[critical_bars,j];
                critical_bars_step(j)=j;
                critical_ten_comp=[critical_ten_comp,2];
                critical_ten_comp_bars(j)=2;
            end
        end
    end
                            
    if isempty(critical_ten_comp)==1
        % if no plastification ocurred in the current loop then the
        % applied loads will be increased by 10 units:
        
        % First the current loads will be removed
        f(edofLoads)=f(edofLoads)-initialLoads;

        initialLoads=initialLoads+dl;
        
        % Then once the current loads are increased, they added to 
        % global force vector
        f(edofLoads)=f(edofLoads)+initialLoads;
        continue;
    else
        collec_dx_nodes=[collec_dx_nodes;
                      Uglobal(1:2:2*nnodes)'];
              
        collec_dy_nodes=[collec_dy_nodes;
                     Uglobal(2:2:2*nnodes)'];
         
        % Computing the unitary deformation of each bar
        
        newCoord(:,1)=coordxy(:,1)+Uglobal(1:2:2*nnodes);
        newCoord(:,2)=coordxy(:,2)+Uglobal(2:2:2*nnodes);

        Lnew=sqrt((newCoord(ni,1)-newCoord(nf,1)).^2+...
             (newCoord(ni,2)-newCoord(nf,2)).^2);
  
        de=(Lnew-L); % Positive when the bar is in tension
        defBarHistory=[defBarHistory,de];
        StrBarHistory=[StrBarHistory,es_bars_normal(1,:)'];
        
        criticalBarsGlobal=[criticalBarsGlobal; 
                              critical_bars_step];
                    
        criticalTenCompBarsGlobal=[criticalTenCompBarsGlobal;
                                      critical_ten_comp_bars];
        
        plasticNsteps=length(StrBarHistory(1,:))-1;
        if plasticNsteps>=2
            for j=1:nbars
                for k=1:plasticNsteps-1
                    if criticalBarsGlobal(k,j)~=0
                        StrBarHistory(j,2+k)=StrBarHistory(j,1+k);
                    end
                end
            end
        end
        
        % Once a bar is plastified by tension or compression it will be
        % substituted by equivalent forces at its nodes corresponding to
        % its resistance force in tension or compression
        A(critical_bars)=1e-15; % its area will be decreased so that the
                                % the stiffness matrix of the whole
                                % structure considers such plastification
        E(critical_bars)=1e-15; % its area will be decreased so that the
                                % the stiffness matrix of the whole
                                % structure considers such plastification
        
        % These are substituiting forces for each plastified bar 
        f(ni(critical_bars)*2-1)=f(ni(critical_bars)*2-1)+...
                                fbe_bars(1,critical_bars);
        f(nf(critical_bars)*2-1)=f(nf(critical_bars)*2-1)+...
                                fbe_bars(3,critical_bars);
        
        f(ni(critical_bars)*2)=f(ni(critical_bars)*2)+...
                            fbe_bars(2,critical_bars);
        f(nf(critical_bars)*2)=f(nf(critical_bars)*2)+...
                            fbe_bars(4,critical_bars);
        
        % Now, once a steel bar has been plastified by tension or
        % compression, it can't suffer another plastification, therefore
        % its resistance is increased to a very high value
        for k=1:length(critical_bars)
            if critical_ten_comp(k)==1
                critical_axial_tension_bar(critical_bars(k))=inf;
            elseif critical_ten_comp(k)==2
                critical_axial_compresion_bar(critical_bars(k))=-inf;
            end
        end
        
        av1=sum(abs(initialLoads1))/length(initialLoads1);
        av2=sum(abs(initialLoads))/length(initialLoads);
        
        FS=av2/av1;
        collection_FS=[collection_FS,FS];
        if iter>1
            % To check the structure's stiffness degradation
            K=zeros(2*nnodes,2*nnodes);
            for i=1:nbars
                [Ke]=bar2e([ex(ni(i)),ex(nf(i))],[ey(ni(i)),ey(nf(i))],...
                    [E(i),A(i)]);
                [K]=assem(Edof(i,:),K,Ke);
            end
            globalKreduced=K(edof,edof);
            det_Kred_post=det(globalKreduced);
            kdam=det_Kred_post/det_Kred;
            collect_Kred=[collect_Kred,det_Kred_post];
            if kdam<kdeg
                break;
            else
                % First the current loads will be removed
                f(edofLoads)=f(edofLoads)-initialLoads;

                initialLoads=initialLoads+dl;

                % Then once the current loads are increased, they added to 
                % global force vector
                f(edofLoads)=f(edofLoads)+initialLoads;
                continue;
            end
        end
    end
end

%%%%%%%%%%%%----------Plotting----------%%%%%%%%%%%%%%%

% Deformation history
if isempty(barsplotPD)==0
    ngt9=0;
    nlet9=0;
    figure(1) % To initialize the plot window
    if barsplotPD(1)>9
        ngt9=ngt9+1;
        bariText2(ngt9,:)=strcat('Bar ',num2str(barsplotPD(1)));
        plot(defBarHistory(barsplotPD(1),:),collection_FS,'k -','LineWidth',1.8)
        legend(bariText2(1,:))
        hold on
    else
        nlet9=nlet9+1;
        bariText(nlet9,:)=strcat('Bar ',num2str(barsplotPD(1)));
        plot(defBarHistory(barsplotPD(1),:),collection_FS,'k -','LineWidth',1.8)
        legend(bariText(1,:))
        hold on
    end
    for i=2:length(barsplotPD)
        if barsplotPD(i)>9
            ngt9=ngt9+1;
            bariText2(ngt9,:)=strcat('Bar ',num2str(barsplotPD(i)));
            
            plot(defBarHistory(barsplotPD(i),:),collection_FS,'LineWidth',1.8,...
                'DisplayName',bariText2(ngt9,:))
        else
            nlet9=nlet9+1;
            bariText(nlet9,:)=strcat('Bar ',num2str(barsplotPD(i)));
            
            plot(defBarHistory(barsplotPD(i),:),collection_FS,'LineWidth',1.8,...
                'DisplayName',bariText(nlet9,:))
        end
    end
    
    xlabel('Delta L')
    ylabel('Qf / Qo')
    title('Deformation history for each bar')
    hold on

    % Stress-Deformation history
    if isempty(barsplotPD)==0
        ngt9=0;
        nlet9=0;
        figure(2)
        if barsplotPD(1)>9
            ngt9=ngt9+1;
        
            bariText2(ngt9,:)=strcat('Bar ',num2str(barsplotPD(1)));
            plot(defBarHistory(barsplotPD(1),:),StrBarHistory(barsplotPD(1),:),...
                'k -','LineWidth',1.8)
            legend(bariText2(1,:))
            hold on
        else
            nlet9=nlet9+1;
            bariText(nlet9,:)=strcat('Bar ',num2str(barsplotPD(1)));
            plot(defBarHistory(barsplotPD(1),:),StrBarHistory(barsplotPD(1),:),...
                'k -','LineWidth',1.8)
            legend(bariText(1,:))
            hold on
        end
        for i=2:length(barsplotPD)
            if barsplotPD(i)>9
                ngt9=ngt9+1;
                bariText2(ngt9,:)=strcat('Bar ',num2str(barsplotPD(i)));

                plot(defBarHistory(barsplotPD(i),:),...
                    StrBarHistory(barsplotPD(i),:),...
                    'LineWidth',1.8,'DisplayName',bariText2(ngt9,:))
            else
                nlet9=nlet9+1;
                bariText(nlet9,:)=strcat('Bar ',num2str(barsplotPD(i)));

                plot(defBarHistory(barsplotPD(i),:),...
                    StrBarHistory(barsplotPD(i),:),...
                    'LineWidth',1.8,'DisplayName',bariText(nlet9,:))
            end
        end

        xlabel('Unitary deformation')
        ylabel('Stress')
        title('Stress-Deformation history for each bar')
        hold on
    end
end

%%------------------- Residual stresses -----------------%%
plasticNsteps=length(StrBarHistory(1,:))-1;
if plasticNsteps>1 && isempty(barsplotPD)==0
    ngt9=0;
    nlet9=0;
    
    residStr(:,1)=zeros(nbars,1);
    residDef(:,1)=zeros(nbars,1);
    residStr(:,2)=StrBarHistory(:,plasticNsteps+1)-StrBarHistory(:,2);
    residDef(:,2)=defBarHistory(:,plasticNsteps+1)-defBarHistory(:,2);
    
    figure(4)
    if barsplotPD(1)>9
        ngt9=ngt9+1;
        bariText2(ngt9,:)=strcat('Bar ',num2str(barsplotPD(1)));
        plot(residDef(1,:),residStr(1,:),'k -','LineWidth',1.8)
        legend(bariText2(1,:))
        hold on
    else
        nlet9=nlet9+1;
        bariText(nlet9,:)=strcat('Bar ',num2str(barsplotPD(1)));
        plot(residDef(1,:),residStr(1,:),'k -','LineWidth',1.8)
        legend(bariText(1,:))
        hold on
    end
    for i=2:length(barsplotPD)
        if barsplotPD(i)>9
            ngt9=ngt9+1;
            bariText2(ngt9,:)=strcat('Bar ',num2str(barsplotPD(i)));
            
            plot(residDef(barsplotPD(i),:),...
                residStr(barsplotPD(i),:),...
                'LineWidth',1.8,'DisplayName',bariText2(ngt9,:))
        else
            nlet9=nlet9+1;
            bariText(nlet9,:)=strcat('Bar ',num2str(barsplotPD(i)));

            plot(residDef(barsplotPD(i),:),...
                residStr(barsplotPD(i),:),...
                'LineWidth',1.8,'DisplayName',bariText(nlet9,:))
        end
    end
    xlabel('Residual deformation')
    ylabel('Residual stress')
    title('Residual stresses-deformation')
    hold on
end

%%------------- Plotting the deformed structure ---------%%

for j=1:nbars
    Ex(j,1)=ex(Edof(j,3)/2);
    Ex(j,2)=ex(Edof(j,5)/2);

    Ey(j,1)=ey(Edof(j,3)/2);
    Ey(j,2)=ey(Edof(j,5)/2);
end

if plotDef==1
    %------------------ Undeformed mesh -----------------%
    figure(3)
    xlabel('Width')
    ylabel('Height')
    title('Deformed-Undeformed Truss');
    plotpar=[1 1 0];
    eldraw2(Ex,Ey,plotpar);

    %---------------------Deformed mesh------------------%
    Ed=extract(Edof,Uglobal); % the last deformation vector before collapse
                              % is considered
    plotpar=[1 3 2];
    eldisp2(Ex,Ey,Ed,plotpar,100);
end

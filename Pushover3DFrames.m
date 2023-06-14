function [barPlasNode,Kfglobal,support,mpbar,plastbars,fglobal]=...
    Pushover3DFrames(qbarxyz,A,Mp,E,G,Iy,Iz,J,coordxyz,ni,nf,eobars,...
    support,bc,seismicforces,dofForces,mpbar,plastbars)

%------------------------------------------------------------------------
% [barPlasNode,Kfglobal,support,mpbar,plastbars,fglobal]=...
% Pushover3DFrames(qbarxyz,A,Mp,E,G,Iy,Iz,J,coordxyz,ni,nf,eobars,...
% support,bc,seismicforces,dofForces,mpbar,plastbars)
%
%------------------------------------------------------------------------
% PURPOSE
%  To compute a static non-linear pushover analysis of a 3D frame
%  
% 
% INPUT:  A = [area_bar;
%               ...]                 area of all elements
%
%         Mp = [Mpi Mpj;             Plastic Moment for each member 
%               ... ]                (i) initial node, (j) final node
%
%         E = [e_bar;                Elasticity modulus of each element
%               ...]                    
%
%         G = [g_bar;                Shear modulus of elasticity of each
%               ...]                 element   
%
%         Iy = [inertia_bar;         momentum of inertia for all elements'
%                       ...]         cross-section with respect to their
%                                    local Y axis (see doc.)
%         
%         Iz = [inertia_bar;         momentum of inertia for all elements'
%                       ...]         cross-section with respect to their
%                                    local Z axis (see doc.)
%
%         J = [polar-inertia_bar;    polar momentum of inertia for all
%                       ...]         elements' cross-sections
%
%         coordxyz = [x,y,z;         node coordinates for all nodes
%                       ...];
%
%         ni                         list of initial nodes of all bars,
%         nf                         list of final nodes of all bars:
%                                         size = [nbars,1]
% 
%         eobars                     local Z axis of each element in the
%                                    global system of reference
%
%         qbarxyz:                   uniformly distributed loads. 
%                                    Size: nbars x 3.
%                                    The first column corresponds to the
%                                    distributed loads in the local X' 
%                                    direction, the second column to the
%                                    loads distributed in the local Y'
%                                    direction and the third column to 
%                                    those in the local Z' direction
%
%         support = [i, j]           support at each bar's end
%                                    options: "Art" or "Fixed"
%                                    (i) initial node, (j) final node
%
%         bc                         restricted dof
%
%         seismicForces = [f(1);]    lateral forces per floor:
%                          f(n);]    size = [nfloors,1]
%
%         Hfloor = [h(1);            Height of each floor from bottom
%                    h(n)]           to top: size = [nfloors,1]
%
%         dofForces = [dof-f(1),     dof at which the lateral forces are
%                       dof-f(n)]    applied (from bottom to top) - global
%
%
%         dydu:                      max ratio between the elastic floor
%                                    deformation and the plastic floor
%                                    deformation -> du/dy -> (dmax-dy)/dy 
%
%         NmaxPlasteps:              Max number of plastic formation steps
%
% OUTPUT: barPlasNode:               vector indicating the plastic
%                                    formations at each element's end at 
%                                    current run:
%                                    1 -> a plastic formation on the
%                                         element's left end
%                                    2 -> a plastic formation on the
%                                         element's right end
%                                    3 -> a plastic formation at each
%                                         element's end
%
%         Kfglobal:                  modified stiffness matrix
%
%         support:                   array indicating the elements' ends 
%                                    conditions
%
%         plastbars:                 array indicating the plastic formation 
%                                    history for each element
%
%         mpbar:                     plastic bending moment caused at the 
%                                    elements' end(s)
%
%         fglobal:                   modified global force vector
%                                    (considering both the seismic inertial
%                                    forces and the internal element
%                                    forces)
% 
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-06-01
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

nbars=length(E);
nnodes=length(coordxyz(:,1));

% Topology matrix 
Edof=zeros(nbars,13);
for i=1:nbars
    Edof(i,1)=i;
    
    Edof(i,2)=ni(i)*6-5;
    Edof(i,3)=ni(i)*6-4;
    Edof(i,4)=ni(i)*6-3;
    Edof(i,5)=ni(i)*6-2;
    Edof(i,6)=ni(i)*6-1;
    Edof(i,7)=ni(i)*6;
    
    Edof(i,8)=nf(i)*6-5;
    Edof(i,9)=nf(i)*6-4;
    Edof(i,10)=nf(i)*6-3;
    Edof(i,11)=nf(i)*6-2;
    Edof(i,12)=nf(i)*6-1;
    Edof(i,13)=nf(i)*6;
    
end

for an=1:2
    Kglobal=zeros(6*nnodes);
    fglobal=zeros(6*nnodes,1);
    fglobal(dofForces)=seismicforces;
        
    elmmat=zeros(12*nbars,12);

    for i=1:nbars      
        
        ex=[coordxyz(ni(i),1) coordxyz(nf(i),1)];
        ey=[coordxyz(ni(i),2) coordxyz(nf(i),2)];
        ez=[coordxyz(ni(i),3) coordxyz(nf(i),3)];
     
        Ex(i,:)=ex;
        Ey(i,:)=ey;
        Ez(i,:)=ez;
        
        eo=eobars(i,:);
        ep=[E(i) G(i) A(i) Iy(i) Iz(i) J(i)];
         
        eq=[0 0 qbarxyz(i,3) 0];

        [Kebar,febar]=beam3e(ex,ey,ez,eo,ep,eq); % This is a CALFEM
                                                 % function
                                                 % Download at: 
                                                 % https://www.byggmek.lth.se/english/calfem/

        if support(i,2)=="Fixed" && support(i,3)=="Art" % DOF: M3
            Mpl=mpbar(i,2);
            eq=[0 0 qbarxyz(i,3) 0 Mpl];
            
            [Kebar,febar]=beamArt3e(ex,ey,ez,ep,eq,eo,1);
             
         elseif support(i,2)=="Art" && support(i,3)=="Fixed" % DOF: M3

             Mpl=mpbar(i,1);
             eq=[0 0 qbarxyz(i,3) 0 Mpl];
             [Kebar,febar]=beamArt3e(ex,ey,ez,ep,eq,eo,2);
         elseif support(i,2)=="Art" && support(i,3)=="Art" % DOF: M3

             Mpl=mpbar(i,1);
             Mp2=mpbar(i,2);
             eq=[0 0 qbarxyz(i,3) 0 [Mpl,Mp2]];
             [Kebar,febar]=beamArt3e(ex,ey,ez,ep,eq,eo,3);
         end
         elmmat((i-1)*12+1:12*i,:)=Kebar; % storing Kebar
         fbe(:,i)=febar; % storing element forces for further use
         
         % Assembling global stiffness matrix
         [Kglobal,fglobal]=assem(Edof(i,:),Kglobal,Kebar,fglobal,febar);

    end 
    
    [Uglobal,Reactions]=solveq(Kglobal,fglobal,bc);
    
    % --- computation of mechanic elements at the ends of bars --- %
    for i=1:nbars
        ue=Uglobal(Edof(i,2:13));
        ke=elmmat((i-1)*12+1:12*i,:);

        fe=-(ke*ue-fbe(:,i));

        reac_bars(:,i)=fe;
    end
    
    if an==1
        current_plas=0; % to register if there is a plastification in the
                        % current load step

        plastified_bars=zeros(nbars,1);
        for i=1:nbars

            % Detect if any end of this bar (i) has been plastified
            bar_plas_check=0;
            if plastbars(1,i)~=0 || plastbars(2,i)~=0
                bar_plas_check=1;
            end
            if bar_plas_check==1
                % Detect if the other end has been also plastified
                if abs(reac_bars(5,i))>=Mp(i,1) && ...
                   abs(reac_bars(11,i))>=Mp(i,2) % if both ends are plastified

                    if plastbars(1,i)==1 && plastbars(2,i)==0 
                        % The bar is currently Art-Fixed and will be Art-Art
                        current_plas=1;

                        plastbars(2,i)=1;

                        % Change condition Fixed-Art to Art-Art
                        support(i,3)="Art";

                        % Storing plastic moment / registration
                        mplas=reac_bars(11,i);
                        mpbar(i,2)=mplas;

                        plastified_bars(i,1)=2;

                    elseif plastbars(1,i)==0 && plastbars(2,i)==1
                        % The bar is currently Fixed-Art and will be Art-Art
                        current_plas=1;
                        plastbars(1,i)=1;

                        % Change condition Fixed-Art to Art-Art
                        support(i,2)="Art";

                        % Storing plastic moment / registration
                        mplas=reac_bars(5,i);
                        mpbar(i,1)=mplas;

                        plastified_bars(i,1)=1;

                    end
                end
            elseif bar_plas_check==0
                if abs(reac_bars(5,i))>=Mp(i,1) && ...
                        abs(reac_bars(11,i))<Mp(i,2)

                    current_plas=1;
                    mplas=reac_bars(5,i);
                    plastbars(1,i)=1;

                    % change condition to Art
                    support(i,2)="Art";

                    % % Storing plastic moment / registration
                    mpbar(i,1)=mplas;

                    plastified_bars(i,1)=1;

                elseif abs(reac_bars(11,i))>=Mp(i,2) && ...
                        abs(reac_bars(5,i))<Mp(i,1)
                    current_plas=1;
                    plastbars(2,i)=1;
                    mplas=reac_bars(11,i);

                    % change condition to Fixed-Art
                    support(i,3)="Art";

                    % Equivalent forces
                    mpbar(i,2)=mplas;

                    plastified_bars(i,1)=2;

                elseif abs(reac_bars(11,i))>=Mp(i,2) && ...
                        abs(reac_bars(5,i))>=Mp(i,1)
                    current_plas=1;
                    plastbars(2,i)=1;
                    plastbars(1,i)=1;

                    mplas1=reac_bars(5,i);
                    mplas2=reac_bars(11,i);

                    % change condition to Art-Art
                    support(i,2)="Art";
                    support(i,3)="Art";

                    % Storing plastic moment / registration
                    mpbar(i,1)=mplas1;
                    mpbar(i,2)=mplas2;

                    plastified_bars(i,1)=3;

                end
            end
        end
    end
    
end
Kfglobal=Kglobal;
barPlasNode=plastified_bars;
%---------------------------------end----------------------------------
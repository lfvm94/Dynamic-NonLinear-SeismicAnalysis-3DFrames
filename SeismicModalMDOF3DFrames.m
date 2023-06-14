function [Cgl,Mgl,Kgl]=SeismicModalMDOF3DFrames(coordxyz,A,unitWeightEl,...
    qbarxyz,eobars,Edof,E,G,J,Iy,Iz,ni,nf,g,AlfaBeta)
% SYNTAX : 
% [Cgl,Mgl,Kgl]=SeismicModalMDOF3DFrames(coordxyz,A,unitWeightEl,...
% qbarxyz,eobars,Edof,E,G,J,Iy,Iz,ni,nf,g,AlfaBeta)
%---------------------------------------------------------------------
%    PURPOSE
%     To compute the global stiffness matrix of a 3D frame as well as
%     the global mass matrix and the global damping matrix with the
%     rayleigh coefficients.
% 
%    INPUT:  coordxyz:          Node coordinates of the structure [x,y,z]
%
%            A:                 Cross-sectional elements' area
%
%            E:                 Modulus of Elasticity of the frame's 
%                               elements
%
%            Iy,Iz:             Cross-sectional inertia with respect to the
%                               local y and z axis of the frame's elements
%
%            bc:                Boundary condition array
%
%            g:                 gravity acceleration
%
%            unitWeightEl:      unit weight material of each element:
%                               Size: nbars x 1
%
%            qbarxyz:           uniformly distributed loads. 
%                               Size: nbars x 3.
%                               The first column corresponds to the
%                               distributed loads in the local X' direction
%                               the second column to the the loads
%                               distributed in the local Y' direction and
%                               the third column to those in the local Z'
%                               direction
%
%    OUTPUT: Cgl:               Global damping matrix
%            Mgl:               Global Mass matrix
%            Kgl:               Global Stiffness matrix
%
%--------------------------------------------------------------------

% LAST MODIFIED: L.Verduzco    2023-06-11
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro
%--------------------------------------------------------------------

nnodes=length(coordxyz(:,1)); nbars=length(E);

%% Stiffness and Mass matrices
Kgl=zeros(6*nnodes);
Mgl=zeros(6*nnodes);
Cgl=zeros(6*nnodes);
for i=1:nbars      
    ex=[coordxyz(ni(i),1) coordxyz(nf(i),1)];
    ey=[coordxyz(ni(i),2) coordxyz(nf(i),2)];
    ez=[coordxyz(ni(i),3) coordxyz(nf(i),3)];
    ep=[E(i) G(i) A(i) Iy(i) Iz(i) J(i)];
    eo=eobars(i,:);
    
    %% Damping Matrix, Stiffness Matrix, Mass matrix
    PV=unitWeightEl(i,1)-qbarxyz(i,3)/A(i); % unit weight of each element
                                          % The distributed downward loads
                                          % on the BEAMS are considered.
    [Mebar,Kebar,Cebar]=FiniteMKCBeams3D(ex,ey,ez,eo,ep,PV,g,AlfaBeta);
    
    [Kgl]=assem(Edof(i,:),Kgl,Kebar);
    [Mgl]=assem(Edof(i,:),Mgl,Mebar);
    [Cgl]=assem(Edof(i,:),Cgl,Cebar);
end 


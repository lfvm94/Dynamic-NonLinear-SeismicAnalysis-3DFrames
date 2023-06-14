% SeismicAnalysis_3DFrames_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To compute the Dynamic Non-Linear Seismic analysis for a 
%    3D Reinforced Concrete Frame.
%
%----------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2023-06-11
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------
clc
clear all

nnodes=8;
nbars=8;

%% Materials
% f'c of each element
fpc=[300;
     300;
     300;
     300;
     300;
     300;
     300;
     300];

% Elasticity modulus of each element in function of f'c
E=14000.*sqrt(fpc);

v=[0.2; % Poisson modulus
    0.2;
    0.2;
    0.2;
    0.2;
    0.2;
    0.2;
    0.2];

G=E./(2.*(1+v)); % Shear modulus

%% Geometry/Topology
% cross-section dimensions of each element (rectangular geometry)
dimensions=[30 30;
            25 40;
            30 30;
            25 40;
            25 40;
            30 30;
            25 40;
            30 30];
        
% cross-section area of each element
A=dimensions(:,1).*dimensions(:,2);
        
% Cross-section inertia
Iy=1/12.*dimensions(:,1).*dimensions(:,2).^3;
Iz=1/12.*dimensions(:,2).*dimensions(:,1).^3;

adim=dimensions(:,2).*0.5;
bdim=dimensions(:,1).*0.5;

% Saint Venant constant (polar inertia - Torsion)
J=adim.*bdim.^3.*(16/3-3.36.*bdim./adim.*(1-bdim.^4./(12.*adim.^4)));
      
%% Topology and node coordinates

% OPTION 1: Manually given
% Coordinates of each node
coordxyz=[0 0 0;
          0 0 300;
          0 500 300;
          0 500 0;
          400 0 0;
          400 0 300;
          400 500 300;
          400 500 0];
      
% Connectivity
NiNf=[1 2;
    2 3;
    3 4;
    3 7;
    2 6;
    5 6;
    6 7;
    7 8];

ni=NiNf(:,1);
nf=NiNf(:,2);

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

%% Prescribed boudnary conditions [dof, displacement]
bc=[1 0;
    2 0;
    3 0;
    4 0;
    5 0;
    6 0;
    19 0;
    20 0;
    21 0;
    22 0;
    23 0;
    24 0;
    25 0;
    26 0;
    27 0;
    28 0;
    29 0;
    30 0;
    43 0;
    44 0;
    45 0;
    46 0;
    47 0;
    48 0];

%% Additional data (optional)
type_elem=[1 "Col";
           2 "Beam";
           3 "Col";
           4 "Beam";
           5 "Beam";
           6 "Col";
           7 "Beam";
           8 "Col"];
       
elemcols=[];
elembeams=[];
beams=0;
cols=0;
for j=1:nbars
    if type_elem(j,2)=="Beam"
        beams=beams+1;
        elembeams=[elembeams,j];
    elseif type_elem(j,2)=="Col"
        cols=cols+1;
        elemcols=[elemcols,j];
    end
end
       
%% Local z axis of each element
eobars=[0 1 0;
        0 0 1;
        0 1 0;
        0 0 1;
        0 0 1;
        0 1 0;
        0 0 1;
        0 1 0];
    
%% Loads       
beams_LL=[1 -40; % Uniformly distributed loads over the beams
          2 -40;
          3 -40;
          4 -40];

% Assignation of distributed loads on beams (local axis [x',y',z'])

qbarxyz=zeros(nbars,4);
qbarxyz(elembeams',3)=beams_LL(:,2);

%% Mode of vibration of interest
modal=[2]; % 2 -> acceleration in the x direction
         % 1 -> acceleration in the y direction
         % [1,2] -> A combination of both (by superposition)
         
%% Damping matrix (for the damped case)
omega1=3;   % Frequencies
omega2=5;

zeta1=0.25;  % Damping factors
zeta2=0.2; 

D=[1/(2*omega1) omega1/2;
    1/(2*omega2) omega2/2];

theta=[zeta1;
    zeta2];

RayCoeff=D\theta; % Rayleigh coefficients

%% Seismic response spectrum from the CFE-15

g=981; % gravity acceleration
Fsit=2.4; FRes=3.8; % Factores de sitio y de respuesta
a0_tau=10; % cm/seg^2

ro=0.8; % Redundance factor
alf=0.9; % Irregularity factor
Q=4; % Seismic behaviour factor

Ta=0.1;
Tb=0.6;
Te=0.5; % Structure's period
k=1.5; % Design spectrum slope
Qp=1+(Q-1)*sqrt(Te/(k*Tb)); % Ductility factor

Ro=2.5; % Over-resistance index
R=Ro+1-sqrt(Te/Ta); % Over-resistance factor

sa=-a0_tau*Fsit*FRes/(R*Qp*alf*ro); % Reduced pseudo-acceleration (cm/seg^2)

%% Modal analysis
pvconc=0.0024; % unit weight of concrete
unitWeightElm=zeros(nbars,1)+pvconc;

% Consistent mass method
[Cgl,Mgl,Kgl]=SeismicModalMDOF3DFrames(coordxyz,A,unitWeightElm,qbarxyz,...
eobars,Edof,E,G,J,Iy,Iz,NiNf(:,1),NiNf(:,2),g,RayCoeff);

Dof=zeros(nnodes,6);
for i=1:nnodes
    Dof(i,1)=6*i-5;
    Dof(i,2)=6*i-4;
    Dof(i,3)=6*i-3;
    Dof(i,4)=6*i-2;
    Dof(i,5)=6*i-1;
    Dof(i,6)=6*i;
end
[Ex,Ey,Ez]=coordxtr(Edof,coordxyz,Dof,2);

%% Dynamic analysis
% Time discretization
dt=0.05;
ttotal=10;
t=0:dt:ttotal;
npoints=length(t);

% Ground acceleration history
tload=1.5; % duration of external excitation

g=sa*cos(5*t); % Acceleration in time
for i=1:length(g)
    if t(i)>tload
        g(i)=20*cos(30*t(i));
    end
end

figure(1)
grid on
plot(t,g,'b -','LineWidth',1.8)
hold on
xlabel('Time (sec)')
ylabel('Acceleration (Kg/cm^2)')
title('Ground acceleration in time')

dofhist=[7]; % dof to evaluate

% Forces history
f(:,1)=zeros(6*nnodes,1);
for i=1:npoints
    % Modal analysis without Damping
    [f(:,i+1),T(:,i),La(:,i),Egv]=ModalsMDOF3DFrames(Mgl,Kgl,...
        bc,g(i),modal);
end

%% Plot of the modal in question and its frequency
Freq=1./T;
zc=0.5*max(coordxyz(:,3));
if length(modal)==1 % If only one modal was entered
    figure(6)
    % Undeformed structure
    grid on
    NoteMode=num2str(modal);
    title(strcat('Eigenmode ','- ',NoteMode))
    elnum=Edof(:,1);
    plotpar=[1,2,1];
    eldraw3(Ex,Ey,Ez,plotpar,elnum)
    
    % Deformed structure
    magnfac=100;
    Edb=extract(Edof,Egv(:,modal));
    plotpar=[1,3,1];
    [magnfac]=eldisp3(Ex,Ey,Ez,Edb,plotpar,magnfac);
    FreqText=num2str(Freq(modal));
    NotaFreq=strcat('Freq(Hz)= ',FreqText);
    text(50,zc,NotaFreq);
end

% Analysis in time with viscous damping
beta=0.25;
gamma=0.5;

d0=zeros(6*nnodes,1);
v0=zeros(6*nnodes,1);

% Solving the motion equation with the Newmark-Beta method
% Initial elements' end support conditions
support=[1 "Fixed" "Fixed";
         2 "Fixed" "Fixed";
         3 "Fixed" "Fixed";
         4 "Fixed" "Fixed";
         5 "Fixed" "Fixed";
         6 "Fixed" "Fixed";
         7 "Fixed" "Fixed";
         8 "Fixed" "Fixed"];
     
% Elastic resistant bending moments for each element's ends
Mp=[50308000 5003800;
    3630000 276940;
    3190000 1190000;
    3000000 57694000;
    4630000 3007694;
    8380000 1380000;
    4630000 4769400;
    9000010 5000010]; % Kg-cm

mpbar=zeros(nbars,2); % to save the plastic moments at each articulation
                      % of each bar as plastifications occur
plastbars=zeros(2,nbars);

% Non-Linear Newmark-Beta
[Dsnap,D,V,A,elPlasHist]=MDOFNewmarkBetaNonLinear3DFrames(Kgl,Cgl,Mgl,d0,...
  v0,dt,beta,gamma,t,f,dofhist,bc,RayCoeff,qbarxyz,A,Mp,E,G,Iy,Iz,J,...
  coordxyz,ni,nf,eobars,support,mpbar,plastbars);

        
%% Dynamic displacement analysis per DOF
figure(2)
grid on
plot(t,D(1,:),'b -','LineWidth',1.8)
legend(strcat('DOF-',num2str(dofhist(1))))
hold on
for i=2:length(dofhist)
    plot(t,D(i,:),'LineWidth',1.8,'DisplayName',...
        strcat('DOF-',num2str(dofhist(i))))
end
hold on
xlabel('Time (sec)')
ylabel('Displacements (cm)')
title('Displacements in time per DOF')

%% Deformation history of structures
dtstep=5;

Xc=max(coordxyz(:,1));
Yc=max(coordxyz(:,2));
Zc=max(coordxyz(:,3));
figure(3)
axis('equal')
axis off
sfac=10;
title(strcat('Deformed structures in time. Scale x ',num2str(sfac)))
for i=1:5
    Ext=Ex+(i-1)*(Xc+400);
    plotpar=[1,2,1];
    elnum=Edof(:,1);
    eldraw3(Ext,Ey,Ez,plotpar,elnum);
    
    Edb=extract(Edof,Dsnap(:,dtstep*i-(dtstep-1)));
    plotpar=[1,3,1];
    eldisp3(Ext,Ey,Ez,Edb,plotpar,sfac);
    Time=num2str(t(5*i-4));
    NotaTime=strcat('Time(seg)= ',Time);
    text((Xc+400)*(i-1)+50,150,NotaTime);
end

Eyt=Ey-(Yc+700);
for i=6:10
    Ext=Ex+(i-6)*(Xc+400);
    plotpar=[1,2,1];
    elnum=Edof(:,1);
    eldraw3(Ext,Eyt,Ez,plotpar,elnum);
    
    Edb=extract(Edof,Dsnap(:,dtstep*i-(dtstep-1)));
    plotpar=[1,3,1];
    [sfac]=eldisp3(Ext,Eyt,Ez,Edb,plotpar,sfac);
    Time=num2str(t(5*i-4));
    NotaTime=strcat('Time(seg)= ',Time);
    text((Xc+400)*(i-6)+50,-250,NotaTime)
    
end
% ----------------------------- End ----------------------------------
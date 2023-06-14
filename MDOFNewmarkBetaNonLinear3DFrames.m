function [Dsnap,D,V,A,elPlasHist]=MDOFNewmarkBetaNonLinear3DFrames(K,C,M,...
    d0,v0,dt,beta,gamma,t,f,dofhist,bc,RayCoeff,qbarxyz,Ae,Mp,Ee,Ge,Iy,Iz,...
    J,coordxyz,ni,nf,eobars,support,mpbar,plastbars)
% SYNTAX : 
% [Dsnap,D,V,A,elPlasHist]=MDOFNewmarkBetaNonLinear3DFrames(K,C,M,...
%  d0,v0,dt,beta,gamma,t,f,dofhist,bc,RayCoeff,qbarxyz,Ae,Mp,Ee,Ge,Iy,Iz,...
%  J,coordxyz,ni,nf,eobars,support,mpbar,plastbars)
%---------------------------------------------------------------------
%    PURPOSE
%     To solve a dynamic system of 2nd order with the numerical
%     integration method "Beta-Newmark".
% 
%    INPUT:  K:                 Global stiffness matrix
%            C:                 Global damping matrix (if any).
%                               Set C=[] for vibration free systems.
%
%            M:                 Global mass matrix
%
%            d0,v0:             Initial displacements and initial
%                               velocities. Vectors of size: n-dof x 1
%
%            beta, gamma:       Chosen parameters for the Beta-Newmark   
%                               method
%
%            t:                 time vector: t0,t1,t2,t3,....tn
%
%            f:                 forces history f(t). Vector of size:
%                               n-dof x n
%
%            dofhist:           degrees of freedom in question
%
%            bc:                boundary condition array. 
%                               Size: n-prescribed-dof x 2
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
%    OUTPUT: D,V,A:             Displacement, velocity and acceleration
%                               history for each DOF in question (taken 
%                               from the vector dofhist)
%
%            Dsnap:             Displacement history for all DOF at each
%                               time step
%
%            elPlasHist:        history of plastic hinge formations at each
%                               element. 
%                               1 -> Plastic formation at the element's
%                               right end
%                               2 -> Plastic formation at the element's
%                               left end
%                               3 -> Plastic formation at both the
%                               element's ends
%
%--------------------------------------------------------------------

% LAST MODIFIED: L.Verduzco    2023-06-13
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro
%--------------------------------------------------------------------

r1=RayCoeff(1); r2=RayCoeff(2); 
[nd,nd]=size(K);
if isempty(C)==1
    C=zeros(nd,nd);  % In case damping is not considered (vibration free)
end
  
nstep=length(t);
b1 = dt^2*0.5*(1-2*beta);  
b2 = (1-gamma)*dt;
b3 = gamma*dt;              
b4 = beta*dt^2;     

ef=(f(:,1)-C*v0-K*d0);
a0=solveq(M,ef,bc); % -> This is a CALFEM function, to download
                    % CALFEM visit its repository at:
                    % https://github.com/CALFEM/calfem-matlab

D(:,1) = d0(dofhist');  % Initial values of solution 
V(:,1) = v0(dofhist');   
A(:,1) = a0(dofhist');       

dnew=d0;    
vnew=v0;    
anew=a0; 
for isnap = 1:nstep
    dpred=dnew+dt*vnew+b1*anew;     
    vpred=vnew+b2*anew;
     
    % Update stiffness matrix
    [elPlasHist(:,isnap),K,support,mpbar,plastbars,f(:,isnap+1)]=...
    Pushover3DFrames(qbarxyz,Ae,Mp,Ee,Ge,Iy,Iz,J,coordxyz,ni,nf,eobars,...
    support,bc,f(:,isnap+1),[1:nd]',mpbar,plastbars);

    % Update Rayleigh Damping Matrix
    C=r1*M+r2*K;
    Keff=M+b3*C+b4*K;
    
    % Update a,v,d for the next iteration
    eff=f(:,isnap+1)-C*vpred-K*dpred;
    anew=solveq(Keff,eff,bc); % -> This is a CALFEM function, to download
                              % CALFEM visit its repository at:
                              % https://github.com/CALFEM/calfem-matlab
    dnew=dpred+b4*anew;  
    vnew=vpred+b3*anew;    

    % Store a,v,d in respective arrays
    D(:,isnap+1) = dnew(dofhist);  % Time-History for each DOF in question
    V(:,isnap+1) = vnew(dofhist);  
    A(:,isnap+1) = anew(dofhist); 

    Dsnap(:,isnap) = dnew; % Time-History displacement for all DOF
    
end
D=D(:,2:nstep+1);
V=V(:,2:nstep+1);
A=A(:,2:nstep+1);
%--------------------------------End----------------------------------
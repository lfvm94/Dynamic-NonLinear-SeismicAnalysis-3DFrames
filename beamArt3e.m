function [Ke,fe]=beamArt3e(ex,ey,ez,ep,eq,eo,art);
% Ke=beamArt3e(ex,ey,ez,ep,eq,eo,art)
% [Ke,fe]=beamArt3e(ex,ey,ez,ep,eq,eo,art,fel)
%----------------------------------------------------------------
%    PURPOSE
%       Calculate the stiffness matrix for a 3D elastic Bernoulli
%       beam element. 
% 
%    INPUT:  ex = [x1 x2]        
%            ey = [y1 y2]   
%            ez = [z1 z2]           node coordinates  
%
%            eo = [xz yz zz];       orientation of local z axis
%
%            ep = [E G A Iy Iz Kv]; element properties 
%                                   E: Young's modulus
%                                   G: Shear modulus 
%                                   A: Cross section area
%                                   Iy: moment of inertia,local y-axis
%                                   Iz: moment of inertia,local z-axis
%                                   Kv: Saint-Venant's torsion constant
% 
%            eq = [qx qy qz qw Mpy];    distributed loads and plastic nodal
%                                       moment with respect to the local Y
%                                       axis
%
%            art = 1, 2 or 3        articulation condition at the ends
%                                   (corresponding to the local 6 and 12 
%                                   DOF) - see documentation:
%                                       1. Fixed - Articulated
%                                       2. Articulated - Fixed
%                                       3. Articulated - Articulated
%
%    OUTPUT: Ke : beam stiffness matrix (12 x 12)
%
%            fe : equivalent nodal forces (12 x 1)
%-----------------------------------------------------------------  

% LAST MODIFIED: L.F.Verduzco    2023-06-01 
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro
%-------------------------------------------------------------


b=[ ex(2)-ex(1); ey(2)-ey(1); ez(2)-ez(1) ];
L=sqrt(b'*b);  n1=b/L;

lc=sqrt(eo*eo'); n3=eo/lc;

qx=eq(1); qy=eq(2); qz=eq(3); qw=eq(4); Mpy=eq(5);
%
E=ep(1); Gs=ep(2);
A=ep(3);
Iy=ep(4); Iz=ep(5);
Kv=ep(6);

a=E*A/L       ; b=12*E*Iz/L^3 ; c=6*E*Iz/L^2;
d=3*E*Iy/L^3 ; e=3*E*Iy/L^2  ; f=Gs*Kv/L;
g=3*E*Iy/L    ; h=2*E*Iz/L    ;
if art==1
    Kle=[a  0  0  0  0  0 -a  0  0  0  0  0 ;
         0  b  0  0  0  c  0 -b  0  0  0  c ;
         0  0  d  0 -e  0  0  0 -d  0  0  0 ;
         0  0  0  f  0  0  0  0  0 -f  0  0 ;
         0  0 -e  0  g 0  0  0  e  0  0  0 ;
         0  c  0  0  0 2*h 0 -c  0  0  0  h ;
        -a  0  0  0  0  0  a  0  0  0  0  0 ;
         0 -b  0  0  0 -c  0  b  0  0  0 -c ;
         0  0 -d  0  e  0  0  0  d  0  0  0 ;
         0  0  0 -f  0  0  0  0  0  f  0  0 ;
         0  0  0  0  0  0  0  0  0  0  0 0 ;
         0  c  0  0  0  h  0 -c  0  0  0 2*h];
     
    fle=L/2*[qx qy 10/8*qz qw qz*L/4 1/6*qy*L qx qy 6/8*qz qw 0 -1/6*qy*L]'+...
        Mpy*[0 0 -3/(2*L) 0 1/2 0 0 0 3/(2*L) 0 1 0]';
elseif art==2
    Kle=[a  0  0  0  0  0 -a  0  0  0  0  0 ;
         0  b  0  0  0  c  0 -b  0  0  0  c ;
         0  0  d  0  0  0  0  0 -d  0 -e  0 ;
         0  0  0  f  0  0  0  0  0 -f  0  0 ;
         0  0  0  0  0  0  0  0  0  0  0  0 ;
         0  c  0  0  0 2*h 0 -c  0  0  0  h ;
        -a  0  0  0  0  0  a  0  0  0  0  0 ;
         0 -b  0  0  0 -c  0  b  0  0  0 -c ;
         0  0 -d  0  0  0  0  0  d  0  e  0 ;
         0  0  0 -f  0  0  0  0  0  f  0  0 ;
         0  0 -e  0  0  0  0  0  e  0  g 0 ;
         0  c  0  0  0  h  0 -c  0  0  0 2*h];

    fle=L/2*[qx qy 6/8*qz qw 0 1/6*qy*L qx qy 10/8*qz qw -qz*L/4 -1/6*qy*L]'+...
        Mpy*[0 0 3/(2*L) 0 1 0 0 0 -3/(2*L) 0 1/2 0]';
elseif art==3
    Mpy1=eq(5);
    Mpy2=eq(6);
    
    Kle=[a  0  0  0  0  0 -a  0  0  0  0  0 ;
         0  b  0  0  0  c  0 -b  0  0  0  c ;
         0  0  0  0  0  0  0  0  0  0  0  0 ;
         0  0  0  f  0  0  0  0  0 -f  0  0 ;
         0  0  0  0  0  0  0  0  0  0  0  0 ;
         0  c  0  0  0 2*h 0 -c  0  0  0  h ;
        -a  0  0  0  0  0  a  0  0  0  0  0 ;
         0 -b  0  0  0 -c  0  b  0  0  0 -c ;
         0  0  0  0  0  0  0  0  0  0  0  0 ;
         0  0  0 -f  0  0  0  0  0  f  0  0 ;
         0  0  0  0  0  0  0  0  0  0  0 0 ;
         0  c  0  0  0  h  0 -c  0  0  0 2*h];

    fle=L/2*[qx qy qz qw 0 1/6*qy*L qx qy qz qw 0 -1/6*qy*L]'+...
        Mpy*[0 0 (Mpy1+Mpy2)/L 0 Mpy1 0 0 0 -(Mpy1+Mpy2)/L 0 Mpy2 0]';
end
%
n2(1)=n3(2)*n1(3)-n3(3)*n1(2);
n2(2)=-n1(3)*n3(1)+n1(1)*n3(3);
n2(3)=n3(1)*n1(2)-n1(1)*n3(2);
%
An=[n1';
    n2;
    n3];
%
G=[  An     zeros(3) zeros(3) zeros(3);
   zeros(3)   An     zeros(3) zeros(3);
   zeros(3) zeros(3)   An     zeros(3);
   zeros(3) zeros(3) zeros(3)   An    ];
%

%
Ke1=G'*Kle*G;  fe1=G'*fle;

Ke=Ke1;
fe=fe1;
%-------------------------- end -------------------------------
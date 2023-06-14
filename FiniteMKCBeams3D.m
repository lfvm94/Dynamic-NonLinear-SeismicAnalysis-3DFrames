function [Me,Ke,Ce]=FiniteMKCBeams3D(ex,ey,ez,eo,ep,PV,ag,AlfaBeta)

% SYNTAX : [Me,Ke,Ce]=FiniteMKCBeams3D(ex,ey,ez,eo,ep,PV,ag,AlfaBeta)
%---------------------------------------------------------------------
%    PURPOSE
%     Compute the mass matrix for a three dimensional beam element. 
% 
%    INPUT:  ex = [x1 x2]
%            ey = [y1 y2]       element node coordinates
%            ez = [z1 z2]
%
%            eo:                local axis of the element in global
%                               coordinates: 
%                               [component-x, component-y, component-z]
% 
%            A                  Transversal area of element
%            PV                 Volumetric weigth of material element
%            ag                 Gravity acceleration

%    OUTPUT: Me : element mass matrix (12 x 12)
%
%--------------------------------------------------------------------
%
% LAST MODIFIED: L.Verduzco    2023-06-10
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro
%--------------------------------------------------------------------

b=[ ex(2)-ex(1); ey(2)-ey(1); ez(2)-ez(1) ];
L=sqrt(b'*b);  n1=b/L;

lc=sqrt(eo*eo'); n3=eo/lc;

%% Stiffness Matrix
E=ep(1); 
Gs=ep(2);
A=ep(3);
Iy=ep(4); Iz=ep(5);
Kv=ep(6);

a=E*A/L       ; b=12*E*Iz/L^3 ; c=6*E*Iz/L^2;
d=12*E*Iy/L^3 ; e=6*E*Iy/L^2  ; f=Gs*Kv/L;
g=2*E*Iy/L    ; h=2*E*Iz/L    ;

Kle=[a  0  0  0  0  0 -a  0  0  0  0  0 ;
     0  b  0  0  0  c  0 -b  0  0  0  c ;
     0  0  d  0 -e  0  0  0 -d  0 -e  0 ;
     0  0  0  f  0  0  0  0  0 -f  0  0 ;
     0  0 -e  0 2*g 0  0  0  e  0  g  0 ;
     0  c  0  0  0 2*h 0 -c  0  0  0  h ;
    -a  0  0  0  0  0  a  0  0  0  0  0 ;
     0 -b  0  0  0 -c  0  b  0  0  0 -c ;
     0  0 -d  0  e  0  0  0  d  0  e  0 ;
     0  0  0 -f  0  0  0  0  0  f  0  0 ;
     0  0 -e  0  g  0  0  0  e  0 2*g 0 ;
     0  c  0  0  0  h  0 -c  0  0  0 2*h];
    
%% Mass Matrix
rx=(Iy+Iz)/A;
R=PV*A*L/(420*ag);
Mle=R*[140  0   0    0     0   0   70   0   0    0     0      0 ;
        0 156   0    0     0 22*L   0  54   0    0     0  -13*L ;
        0   0   156  0  -22*L  0    0   0   54   0  13*L      0 ;
        0   0   0 140*rx^2 0   0    0   0   0 70*rx^2  0      0 ;
        0   0 -22*L  0  4*L^2  0    0   0 -13*L  0 -3*L^2     0 ;
        0 22*L  0    0     0 4*L^2  0 13*L  0    0     0 -3*L^2 ;
       70   0   0    0     0   0   140  0   0    0     0      0 ;
        0  54   0    0     0 13*L   0 156   0    0     0  -22*L ;
        0   0   54   0  -13*L  0    0   0  156   0  22*L      0 ;
        0   0   0 70*rx^2  0   0    0   0   0 140*rx^2 0      0 ;
        0   0  13*L  0  -3*L^2 0    0   0 22*L   0  4*L^2     0 ;
        0 -13*L 0    0     0 -3*L^2 0 -22*L 0    0     0  4*L^2];


%% Damping Matrix
a=AlfaBeta(1); b=AlfaBeta(2);
Cle=a*Mle+b*Kle;

%% Rotation to global coordinates
n2(1)=n3(2)*n1(3)-n3(3)*n1(2);
n2(2)=-n1(3)*n3(1)+n1(1)*n3(3);
n2(3)=n3(1)*n1(2)-n1(1)*n3(2);

An=[n1';
    n2;
    n3];

G=[  An     zeros(3) zeros(3) zeros(3);
   zeros(3)   An     zeros(3) zeros(3);
   zeros(3) zeros(3)   An     zeros(3);
   zeros(3) zeros(3) zeros(3)   An    ];
  
 Me=G'*Mle*G;   Ke=G'*Kle*G;  Ce=G'*Cle*G;
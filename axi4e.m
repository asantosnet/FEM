function [Ke,fe]=axi4e(ex,ey,ep,D,eq)
%Ke=axi4e(ex,ey,ep,D)
%[Ke,fe]=axi4e(ex,ey,ep,D,eq)
%
% Purpose: Compute element stiffness matrix (and load vector)
%          for a 4-node quadrilateral isoparametric element
%          The z-axis corresponds to a circumferential direction  
% Input:   ex = [x1 x2 x3 x4] - x-coordinates for the element
%          ey = [y1 y2 y3 y4] - y-coordinates for the element
%          ep = [xa ya nax nay ngauss]
%                xa     - x-coordinate on the axis of rotation
%                ya     - y-coordinate on the axis of rotation
%                nax    - x-component of a direction vector thru (xa,ya)
%                         describing the axis orientation
%                nay    - y-component of a direction vector thru (xa,ya)
%                         describing the axis orientation
%                ngauss - number of gauss points in each direction
%          D  - constitutive matrix for the axi-symmetric case
%          eq = [bx by]'
%                bx     - x-component of bodyforce 
%                by     - y-component of bodyforce

% Element routine developed September 2006  
%   Fredrik Larsson, 
%   Chalmers University of Technology,
%   fredrik.larsson@chalmers.se
% Revisions
%-----------
% R1    October 18, 2006    Error in intregrationscheme corrected
%
if nargin<5
    eq=[0;0];
end

nz=[ep(3) ep(4)]'./sqrt(ep(3).^2+ep(4).^2);
nr=[nz(2) -nz(1)]';
L=[nr nz]';

erz=L*([ex-ep(1);ey-ep(2)]);
er=erz(1,:);
ez=erz(2,:);

Leff=zeros(8,8);
for i=1:4
    Leff(2*(i-1)+[1 2],2*(i-1)+[1 2])=L;
end

ngauss=ep(5);

if ngauss==1
    gauss=[1;2];
elseif ngauss==2
    gauss=[-.577350269189626 .577350269189626;1 1];
else
    if ngauss~=3
        disp('Requested number of gausspoints not implemented.')
        disp('Uses 3x3 gausspoints.')
    end
    gauss=[-0.774596669241483 0 0.774596669241483;(2-.888888888888889)/2 .888888888888889 (2-.888888888888889)/2];
end

Ke=zeros(8,8);
fe=zeros(8,1);

for i=1:ngauss
    for j=1:ngauss
        xi=gauss(1,i);
        eta=gauss(1,j);
        H=gauss(2,i).*gauss(2,j);
        N=[(xi-1)*(eta-1) -(xi+1)*(eta-1) (xi+1)*(eta+1) -(xi-1)*(eta+1)]./4;
        Ne=zeros(2,8);
        Ne(1,1:2:7)=N;
        Ne(2,2:2:8)=N;
        
        r=N*er';
        
        GradN=[eta-1 -(eta-1) eta+1 -(eta+1);...
                      xi-1  -(xi+1)  xi+1  -(xi-1)];
        J=[GradN(1,:)*er' GradN(2,:)*er';...
            GradN(1,:)*ez' GradN(2,:)*ez'];
        GradNxy=J'\GradN;
        
        Be=zeros(4,8);
        Be(1,1:2:7)=GradNxy(1,:);
        Be(2,2:2:8)=GradNxy(2,:);
        Be(3,1:2:7)=N./r;
        Be(4,1:2:7)=GradNxy(2,:);
        Be(4,2:2:8)=GradNxy(1,:);        

        Ke=Ke+Be'*D*Be.*r.*H.*det(J);
        fe=fe+Ne'*eq.*r.*H.*det(J);
    end    
end

Ke=Leff'*Ke*Leff;
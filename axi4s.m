function [es,et,eci]=axi4s(ex,ey,ep,D,ed)
%[es,et,eci]=axi4s(ex,ey,ep,D,ed)
%
% Purpose: Compute stresses and strains in an axisymmetric 
%          4-node quadrilateral isoparametric element
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
%          ed = [ux1 uy1 ux2 uy2 ux3 uy3 ux4 uy4] - element displacements
% Output:  es =  [sigma^1_xx sigma^1_yy sigma^1_zz sigma^1_xy;
%                 sigma^2_xx sigma^2_yy sigma^2_zz sigma^2_xy;
%                 ..         ..         ..         ..        ]
%                - stresses where each row corresponds du one gausspoint
%          et  = [epsilon^1_xx epsilon^1_yy epsilon^1_zz gamma^1_xy;
%                 epsilon^2_xx epsilon^2_yy epsilon^2_zz gamma^2_xy;
%                 ..         ..         ..         ..              ]
%                - strains where each row corresponds du one gausspoint
%          eci = [x_1 y_1;
%                 x_2 y_2;
%                 ..  .. ]
%                - coordinates for each gausspoints

% Element routine developed September 2006  
%   Fredrik Larsson, 
%   Chalmers University of Technology,
%   fredrik.larsson@chalmers.se

nz=[ep(3) ep(4)]'./sqrt(ep(3).^2+ep(4).^2);
nr=[nz(2) -nz(1)]';
L=[nr nz]';

erz=L*([ex-ep(1);ey-ep(2)]);
er=erz(1,:);
ez=erz(2,:);

edrz=L*[ed(1:2:7);ed(2:2:8)];
edrznew=zeros(1,8);
edrznew(1:2:7)=edrz(1,:);
edrznew(2:2:8)=edrz(2,:);

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

es=zeros(ngauss^2,4);
et=zeros(ngauss^2,4);
eci=zeros(ngauss^2,2);

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
        z=N*ez';
        
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

        ig=(i-1)*ngauss+j;
        
        epsrz=Be*edrznew';
        eps=L*[epsrz(1) epsrz(4);epsrz(4) epsrz(2)]*L';
        et(ig,:)=[eps(1,1) eps(2,2) epsrz(3) eps(1,2)];
        
        sigrz=D*epsrz;
        sig=L*[sigrz(1) sigrz(4); sigrz(4) sigrz(2)]*L';
        es(ig,:)=[sig(1,1) sig(2,2) sigrz(3) sig(1,2)];
        
        eci(ig,:)=([ep(1);ep(2)]+L'*[r;z])';
    end    
end
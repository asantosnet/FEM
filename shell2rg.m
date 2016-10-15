function [ KMe,fe,KGe] = shell2rg( ex,ey,ep,D,es,G,plate,eq)
%Purpose : Compute the Geometrical stiffness matrix, global material matrix
%and the force fe matrix for Blucking
%
% Input
% ex = [x1 x2 x3 x4]     element coordinates
% ey = [y1 y2 y3 y4]
% ep = [ptype t ]             ptype: analysis type,    t: thickness
% G 
% D                              constitutive matrix
% ed = [u1 u2 .. u8;          element displacement vector
%                ...........]          one row for each element
% plate                     Minldlin if 1 and Kirchoff if 2
% eq =[qx qy qz]
%
% Output
% KMe = 20x20               
% fe = 20x1                 Global material matrix
% KGe = 20x20               Geometrical stiffness matrix,
% ec = [Kxx Kyy 2Kxy]     curvature

sigma = [ es(1) es(3);
    es(3) es(2)];

% C matrix

a = (1/2)*(ex(3)-ex(1));
b = (1/2)*(ey(3)-ey(1));

C = [1 -a -b  a^2  a*b  b^2     -a^3 -(a^2)*b -(b^2)*a     -b^3     (a^3)*b    a*(b^3);
    0  0  1    0   -a  -2*b        0     (a^2)   2*b*a   3*(b^2)     -(a^3) -3*a*(b^2);
    0 -1  0  2*a    b     0  -3*(a^2)  -2*a*b    -(b^2)       0   3*(a^2)*b      (b^3);
    1  a -b  a^2 -a*b   b^2      a^3 -(a^2)*b  (b^2)*a     -b^3    -(a^3)*b   -a*(b^3);
    0  0  1    0    a  -2*b        0     (a^2)  -2*b*a  3*(b^2)       (a^3)  3*a*(b^2);
    0 -1  0 -2*a    b     0 -3*(a^2)    2*a*b    -(b^2)       0   3*(a^2)*b      (b^3);
    1  a  b  a^2  a*b   b^2      a^3  (a^2)*b  (b^2)*a      b^3     (a^3)*b    a*(b^3);
    0  0  1    0    a   2*b        0     (a^2)   2*b*a  3*(b^2)       (a^3)  3*a*(b^2);
    0 -1  0 -2*a   -b     0 -3*(a^2)   -2*a*b    -(b^2)       0  -3*(a^2)*b     -(b^3);
    1 -a  b  a^2 -a*b   b^2     -a^3  (a^2)*b -(b^2)*a      b^3    -(a^3)*b   -a*(b^3);
    0  0  1    0   -a   2*b        0     (a^2)  -2*b*a  3*(b^2)      -(a^3) -3*a*(b^2);
    0 -1  0  2*a   -b     0 -3*(a^2)    2*a*b    -(b^2)       0  -3*(a^2)*b     -(b^3)];

% Gauss points for a 3x3 integration

xGauss = [0 0.774596669241483  -0.774596669241483]*a;
yGauss = [0 0.774596669241483 -0.774596669241483]*b;

Hx = [0.888888888888889 0.555555555555556   0.555555555555556]*a;
Hy = [0.888888888888889 0.555555555555556   0.555555555555556]*b;

% calculating deltaN , SEE IF THE -1 IS CORRECT IN THIS WAY

KGeep = zeros(12,12);

% Calculatig the components of KGeep for each gauss points in a 3x3
% integration

for i = 1:3
    for j = 1:3
        
        % Calculating the delta*N
        
        deltaN1 = [0 1 0 2*xGauss(i) yGauss(j) 0 3*(xGauss(i)^2) 2*xGauss(i)*yGauss(j) yGauss(j)^2 0 3*(xGauss(i)^2)*yGauss(j) yGauss(j)^3];
        
        deltaN2 = [0 0 1 0 xGauss(i) 2*yGauss(j) 0 (xGauss(i)^2) 2*yGauss(j)*xGauss(i) 3*(yGauss(j)^2) (xGauss(i)^3) 3*xGauss(i)*(yGauss(j)^2)];
        
        deltaN = [deltaN1;deltaN2];
        
        
        deltaN = deltaN*(C^-1);
        
        Nr = ep(2)*sigma;
        
        KGeep = KGeep + deltaN'*Nr*deltaN*Hx(i)*Hy(j);
        
    end
    
end

posep = [3 4 5 8 9 10 13 14 15 18 19 20];

% Putting KGeep in their proper palces
KGe(posep,posep) = KGeep;

% Calculating KMe and fe
[KMe,fe] = shell2re(ex,ey,ep,D,G,eq,plate);

end

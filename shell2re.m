function [ Ke,Fe ] = shell2re(ex,ey,ep,D,G,eq,plate)
% 
% Purpose : Calculate the Element stiffness matrix and the Element load
% vector
% 
% --------------------------------------------------------
% OUTPUT -
%
% Ke = Element stiffness matrix (20 x 20)
% Fe = ELement load vector (20 x 1)
%
% INPUT -
% ex = [x1 x2 x3 x4] X -Coordinates
% ey = [x1 x2 x3 x4] Y- coordinates
% ep = [ptype,t,irb,irs] - ptype = 1 gives plane stress
%                                  2 gives plane strain
%                          t thickness of the plate
%                          irb and irs integration rules for bending and
%                          shearing respectively (only for Mindlin plate)
% D = Constitutive matrix
% G = Constitutive matrix for bending mode (Mindlin plate) 
% eq = [bx,by,qz]
% plate = 1 for Kirchoff plate and 2 for Mindlin plate

Ke = zeros (20,20);
Fe = zeros (20,1);
[keplanre,feplanre] = planre ([ex(1) ex(3)],[ey(1) ey(3)],[ep(1) ep(2)],D,[eq(1) eq(2)]);

% TO swicth betwehn Krichof and Mildlin theory
switch plate
    
    case 1 % Kirchoff theory
        [kep,fep] = platre (ex,ey,ep(2),D,eq(3));

    case 2 % Mindlin theory        
        [kep,fep] = platme (ex,ey,[ep(2) ep(3) ep(4)],D,G,eq(3));      
end

% To add Fe and Ke the good position - Easier way exits

posplanre = [1 2 6 7 11 12 16 17];
posep = [3 4 5 8 9 10 13 14 15 18 19 20];

Ke(posep,posep) = kep; 
Ke(posplanre,posplanre) = keplanre;
Fe(posep) = fep ;
Fe(posplanre) =feplanre;

end


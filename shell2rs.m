function [es,et,ef,ec] = shell2rs(ex,ey,ep,D,ed)

% Purpose : Compute element stresses, strains, forces and curvature
%
% Input
% ex = [x1 x2 x3 x4]     element coordinates
% ey = [y1 y2 y3 y4]
% ep = [ptype t ]             ptype: analysis type,    t: thickness                                 
% D                           constitutive matrix
% ed = [u1 u2 .. u8;          element displacement vector
%                ...........]          one row for each element
% 
% Output
% es = [sigxx sigyy sigxy]    stresses
% et = [epsxx epsyy gamxy]     strains
% ef = [ Mxx Myy Mxy Vxz Vyz] section forces
% ec = [Kxx Kyy 2Kxy]     curvature


% element force matrix and curvature with platrs
[ef,ec]=platrs(ex,ey,ep(2),D,ed([3 4 5 8 9 10 13 14 15 18 19 20]));

% element stresses and strains with planrs
[es,et]=planrs([ex(1) ex(3)],[ey(1) ey(3)],[ep(1) ep(2)],D,ed([1 2 6 7 11 12 16 17]));



end


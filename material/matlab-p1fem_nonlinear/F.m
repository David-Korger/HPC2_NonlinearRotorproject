function Fu = F(u,param)
%*** Assembly of stiffness matrix
d = zeros(size(param.coordinates,1),1);
for j = 1:size(param.elements,1)
  nodes = param.elements(j,:);
  vertices = param.coordinates(nodes,:);
  G = [ones(1,3);vertices'] \ [zeros(1,2);eye(2)];
  Du = G'*u(nodes);
  xM = sum(vertices)/3;
  cA = OpA(xM,Du,param.material(j));
  d(nodes) = d(nodes) ...
  + 1/2*det([ones(1,3);vertices']) * G * [cA(1:2);cA(2:3)] * Du;
end
%*** Computation of F(u) = d(u) - b  
Fu = zeros(size(param.coordinates,1),1);
Fu(param.freenodes) = d(param.freenodes) - param.b(param.freenodes);



% %*** Assembly of stiffness matrix
% A = sparse(size(param.coordinates,1),size(param.coordinates,1));
% for j = 1:size(param.elements,1)
%   vertices = param.coordinates(param.elements(j,:),:);
%   G = [ones(1,3);vertices'] \ [zeros(1,2);eye(2)];
%   xM = sum(vertices)/3;
%   gradu = G'*u(param.elements(j,:));
%   cA = OpA(xM,gradu,param.material(j));
%   A(param.elements(j,:),param.elements(j,:)) = ...
%     A(param.elements(j,:),param.elements(j,:)) ...
%   + 1/2*det([ones(1,3);vertices']) * G * [cA(1:2);cA(2:3)] * G';
% end
% %*** Computation of F(u) = b -A(u) * u 
% Fu = zeros(size(param.coordinates,1),1);
% Fu(param.freenodes) = A(param.freenodes,:) * u - param.b(param.freenodes);

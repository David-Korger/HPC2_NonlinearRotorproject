function DFu = DF(u,para)
%*** Assembly of Jacobian matrix
DFu = sparse(size(para.coordinates,1),size(para.coordinates,1));
for j = 1:size(para.elements,1)
  vertices = para.coordinates(para.elements(j,:),:);
  G = [ones(1,3);vertices'] \ [zeros(1,2);eye(2)];
  [val,Dval] = OpA(sum(vertices)/3,G'*u(para.elements(j,:)),para.material(j));
  DFu(para.elements(j,:),para.elements(j,:)) = ...
    DFu(para.elements(j,:),para.elements(j,:))  ...
  + 1/2*det([ones(1,3);vertices']) * G * [Dval(1,1:2);Dval(1,2:3)] * G';
end

% for j = 1:size(para.elements,1)
%   vertices = para.coordinates(para.elements(j,:),:);
%   G = [ones(1,3);vertices'] \ [zeros(1,2);eye(2)];
%   Du = G'*u(para.elements(j,:));
%   absDu = norm(Du);
%   M = p(absDu,para.material(j)) * eye(2);
%   if absDu >1e-10; 
%     M = M + Dp(absDu,para.material(j))/absDu * (Du * Du');
%   end
%   DFu(para.elements(j,:),para.elements(j,:)) = ...
%     DFu(para.elements(j,:),para.elements(j,:))  ...
%   + 1/2*det([ones(1,3);vertices']) * G * M * G';
% end

function A = DFsimple(u,param)
A = sparse(size(u,1),size(u,1));
delta = 1e-5;
ej = zeros(size(u,1),1);
for j=1:size(u,1)
  ej(j) = 1;
  A(:,j) = (F(u+delta*ej,param) - F(u-delta*ej,param))/(2*delta); 
  ej(j) = 0;
end
clear all
% FEM2D  two-dimensional finite element method for linear second order PDE.
% Initialisation
load coordinates.dat; 
load elements.dat;
eval('load neumann.dat;','neumann=[];');
eval('load material.dat;','material=ones(size(elements,1),1);');
load dirichlet.dat; 
%*** Graphic representation
figure(1)
patch(reshape(coordinates(elements,1),[],3)', ...
      reshape(coordinates(elements,2),[],3)',material(:,[1,1,1])')
view(2)
%*** Refinement
for k=1:4
  [coordinates,elements,material,dirichlet,neumann] ...
     = refineR(coordinates,elements,material,dirichlet,neumann);
end
%*** Prescribe values at Dirichlet nodes
para.dirichlet = unique(dirichlet);
u = zeros(size(coordinates,1),1);
u(para.dirichlet) = uD(coordinates(para.dirichlet,:)); 
%*** Store mesh as structure
para.elements = elements;
para.coordinates = coordinates;
para.material = material;
para.freenodes = setdiff(1:size(coordinates,1), para.dirichlet);
%*** Right hand side
para.b = zeros(size(coordinates,1),1);
for j = 1:size(elements,1)
  vertices = coordinates(elements(j,:),:);  
  xM = sum(vertices)/3;
  para.b(elements(j,:)) = para.b(elements(j,:)) ...
   +det([1 1 1;vertices']) * f1(xM,material(j))/6 ...
   -1/2*(vertices([2,3,1],:)-vertices([3,1,2],:))*[0,-1;1,0]*f2(xM,material(j))';
end
%*** Neumann conditions
for j = 1 : size(neumann,1)
   n = coordinates(neumann(j,2),:) - coordinates(neumann(j,1),:);
   dist = norm(n);
   n = ([0,1;-1,0]*n')/dist;
   xM = sum(coordinates(neumann(j,:),:))/2;
   para.b(neumann(j,:)) = para.b(neumann(j,:)) ...
                        + 0.5 * dist * (g(xM)+f2(xM,material(j))*n);
end
%*** Solve with Newton and variable step size
[u,nit] = newton(u,@F,@DF,1e-10,20,para);
fprintf('No. of Newton iterations = %3i\n',nit)
%*** Graphic representation
figure(2)
trisurf(elements,coordinates(:,1),coordinates(:,2),u,'facecolor','interp')

view(-16,20)
print -depsc2 'Solution_nonlinear'
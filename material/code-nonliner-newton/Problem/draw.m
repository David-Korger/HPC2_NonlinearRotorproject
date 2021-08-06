function draw
clf
C = load('motor.co');
nC = size(C,1);

% data = load('motor.bd');
% D = [data(:,1),data(:,3)+nC;data(:,3)+nC,data(:,2)]+1;

data = load('motor.el');
nE = max(max(data(:,4:6))+1);
C = [C; zeros(nE,2)];
C(nC+data(:,4)+1,:) = (C(data(:,1)+1,:)+C(data(:,2)+1,:))/2;
C(nC+data(:,5)+1,:) = (C(data(:,2)+1,:)+C(data(:,3)+1,:))/2;
C(nC+data(:,6)+1,:) = (C(data(:,3)+1,:)+C(data(:,1)+1,:))/2;

E = [data(:,1),data(:,4)+nC,data(:,6)+nC; ...
     data(:,2),data(:,5)+nC,data(:,4)+nC; ...
     data(:,3),data(:,6)+nC,data(:,5)+nC; ...
     data(:,5)+nC,data(:,6)+nC,data(:,4)+nC]+1;
Material = repmat(data(:,7),4,1);
nC,nE

x = load('motor.sol');
u = x(1:size(C,1));

figure(1)
trisurf(E,C(:,1),C(:,2),u,'facecolor','interp')
% hold on
% plot(reshape(C(D,1),[],2), reshape(C(D,2),[],2), '*-');
% hold off
axis equal; axis off, view(2)

M = 230;
N = 40;
d = 0.1;

s = linspace(-d+min(C(:,1)),d+max(C(:,1)),M);
t = linspace(-d+min(C(:,2)),d+max(C(:,2)),M);

[x,y] = meshgrid(s,t);
z = tri2monic(C,E,u,x,y);
z = (z{1}+min(z{1}(:)))./(max(z{1}(:))-min(z{1}(:)));
figure(2)
surf(x,y,z);
figure(3)
contour(x,y,z,N,'k')
hold on
trisurf(E,C(:,1),C(:,2),-0.1+0*u,Material,'edgecolor','none')
hold off
view(2)
axis equal
axis off




function ugrid = tri2monic(coordinates,elements,u,x,y)
% interpolate piecewise affine function onto grid (x,y)

sx = x(1,:); sy = y(:,1);
m = size(x,2); n = size(x,1);
hx = sx(2)-sx(1); hy = sy(2)-sy(1);

ugrid = cell(size(x,1), size(x,2), size(u,2));
for j = 1:size(u,2)
  ugrid{j} = NaN(size(x,1),size(x,2));
end

% Gehe alle Dreiecke durch
for j = 1:size(elements,1)
     p = coordinates(elements(j,:),:);
        
     % Gitter in X-Richtung bestimmen
     xminDreieck = max(1,floor((min(p(:,1)) - sx(1))/hx)+1);
     xmaxDreieck = min(m,ceil( (max(p(:,1)) - sx(1))/hx));
     dsx = sx(xminDreieck:xmaxDreieck);
     
     % Gitter in Y-Richtung bestimmen
     yminDreieck = max(1,floor((min(p(:,2)) - sy(1))/hy)+1);
     ymaxDreieck = min(n,ceil( (max(p(:,2)) - sy(1))/hy));
     dsy = sy(yminDreieck:ymaxDreieck);
     
     % Gitter erzeugen
     [dx,dy] = meshgrid(dsx,dsy);
     dugrid = cell(1, size(u(elements(j,:)),2));
     for k = 1:size(u(elements(j,:)),2)
         dugrid{k} = NaN(size(dx,1),size(dx,2));
     end
     % tri2grid aufrufen
     dugrid = tri2grid(p,u(elements(j,:)),dx,dy,dugrid);
     % das kleine Grid in das groessere einfuegen
     for k = 1:size(u(elements(j,:)),2)
         ugrid{k}(yminDreieck:ymaxDreieck,xminDreieck:xmaxDreieck) = ...
             min(dugrid{k},ugrid{k}(yminDreieck:ymaxDreieck,xminDreieck:xmaxDreieck));
     end
end


function ugrid = tri2grid(P,u,x,y,ugrid)
% interpolate from triangle with vertices P to monic grid (x,y)
% x, y, and all ugrid's{:} have to have the same size

[m,n] = size(x);
x = x(:); y = y(:);
L = zeros(m*n,3);
L(:,1) = (P(2,1)-x).*(P(3,2)-y) - (P(3,1)-x).*(P(2,2)-y);
L(:,2) = (P(3,1)-x).*(P(1,2)-y) - (P(1,1)-x).*(P(3,2)-y);
L(:,3) = (P(1,1)-x).*(P(2,2)-y) - (P(2,1)-x).*(P(1,2)-y);

ind = find (L(:,1)>=0 & L(:,2)>=0 & L(:,3)>=0);

if ~isempty(ind)
  area = sum(L(ind,:),2);
  for j = 1:size(u,2)
    ugrid{j}(ind) = (L(ind,:)*u)./area;
    ugrid{j} = reshape(ugrid{j},m,n);
  end
end




function [coordinates,newElements,newMaterial,varargout] ...
           = refineR(coordinates,elements,material,varargin)
%*** Obtain geometric information on edges
%    Collect all edges
edges = reshape(elements(:,[1,2,3,2,3,1]),[],2);          
ptr = [3*size(elements,1),zeros(1,nargin-3)]; 
for j = 1:nargin-3
  ptr(j+1) = size(varargin{j},1);
  edges = [edges;varargin{j}];
end   
ptr = cumsum(ptr);
%*** Create numbering of edges
[edge2nodes,~,ie] = unique(sort(edges,2),'rows');
element2edges = reshape(ie(1:ptr(1)),[],3);
%*** Provide boundary2edges
for j = 1:nargin-3
  boundary2edges{j} = ie(ptr(j)+1:ptr(j+1));
end
%*** Generate new nodes
edge2newNode = size(coordinates,1) + (1:max(element2edges(:)))';
coordinates = [coordinates; (coordinates(edge2nodes(:,1),:) ...
                           + coordinates(edge2nodes(:,2),:))/2];
%*** Refine boundary conditions
for j = 1:nargout-3
  varargout{j} = [varargin{j}(:,1),edge2newNode(boundary2edges{j}); ...
                  edge2newNode(boundary2edges{j}),varargin{j}(:,2)];
end
%*** Provide new nodes for refinement of elements
nodes = [elements,edge2newNode(element2edges)];
%*** Generate new material
material = material(:);
newMaterial = reshape(material(:,[1,1,1,1]),[],1);
%*** Generate new elements
newElements = [nodes(:,[1,4,6]);nodes(:,[4,2,5]);nodes(:,[6,5,3]);nodes(:,[5,6,4])];
     
function M = stima(vertices,material)
d = size(vertices,2);
G = [ones(1,d+1);vertices'] \ [zeros(1,d);eye(d)];
xM = sum(vertices)/(d+1);
A = OpA(xM,[],material);
M = det([ones(1,d+1);vertices']) * G * [A(1:2);A(2:3)] * G' / prod(1:d);
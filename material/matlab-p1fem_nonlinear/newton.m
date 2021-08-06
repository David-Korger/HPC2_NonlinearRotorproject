function [u,nit] = newton(u,F,DF,tol,maxit,param)
Fu = F(u,param);
DFu = DF(u,param);
u',DFu
s = zeros(size(u,1),1);
s(param.freenodes) = -DFu(param.freenodes,param.freenodes)\Fu(param.freenodes);
lam = 1;
tmp = max(tol,tol*norm(s));
nit=0;
while norm(s) > tmp && nit <= maxit
  nit = nit + 1;
  u_old = u;
  lam = min(1,2*lam);
  for k=1:30
    u = u_old + lam * s; % Daempfung mit Parameter lam
    Fu = F(u,param);
    s0 = -DFu(param.freenodes,param.freenodes)\Fu(param.freenodes);
    if norm(s0) <= (1-lam/2) * norm(s) % Abbruch der for-Schleife, falls
      break                            % Konvergenztest erfuellt
    end                                % lam noch zu groß --> halbieren
    lam = lam/2;
  end
  DFu = DF(u,param);
  s(param.freenodes) = -DFu(param.freenodes,param.freenodes)\Fu(param.freenodes);
end


function [val,Dval] = OpA(x,Du,material)
absDu = norm(Du);
val = [p(absDu,material),zeros(size(x,1),1),p(absDu,material)];
if nargout>1
  Dval = val;
  if absDu > 1e-12
    Dval = Dval + Dp(absDu,material)*...
                 [Du(1)^2,Du(1)*Du(2),Du(2)^2];
  end
end

function value = p(t,material)
if material ==1
  value = 1;
elseif material == 2
  value = 2-1/(1+t);
end

function value = Dp(t,material)
if material ==1
  value = 0;
elseif material == 2
  value = +1./(t*(1+t).^2);
end
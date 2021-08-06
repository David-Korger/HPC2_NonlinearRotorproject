function value = p(t,material)
if material ==1
  value = 2+1./(1+t(:));
elseif material == 2
  value = 4+1./(1+t(:));
else
  value = .5+1./(1+t(:));
end
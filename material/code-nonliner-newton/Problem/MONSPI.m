clear all, format compact, format long g, clf
H = [ 400,  600,  800, 1000, 1400, 2000, 3000, 4000, 6000];
B = [0.73, 0.92, 1.05, 1.15, 1.28, 1.42, 1.53, 1.565, 1.60];
nu_inf = 24000;

H = [20 60 80 95 105 120 140 160 180 200 240 2500];
B = [0.19 0.65 0.87 1.04 1.18 1.24 1.272 1.3 1.32 1.34 1.36 1.45];

nu = H./B;
nu_inf = (H(end)-H(end-1))/(B(end)-B(end-1));
% *** Step 1, get slope on each element and average on neighbours
Delta = (nu(2:end)-nu(1:end-1))./(B(2:end)-B(1:end-1)); 
d =  [0,( (B(2:end-1)-B(1:end-2)).*Delta(2:end) ...
        + (B(3:end  )-B(2:end-1)).*Delta(1:end-1) )./(B(3:end)-B(1:end-2))];
d = [d, (- (nu(end)-nu(end-2)) * (B(end)-B(end-1))./(B(end)-B(end-2)) ...
        + Delta(end)* (B(end)-B(end-2)))/(B(end-1)-B(end-2))];
% *** Step 2, reduce slopes
d(end) = min(max(0,d(end)),2*(nu_inf-nu(end))/B(end));
a = d(1:end-1)./Delta;
b = d(2:end)./Delta;
c = 1;
for k=1:length(Delta)
  d(k) = d(k) * c;
  a(k) = a(k) * c;
  c = min(3/sqrt(a(k)^2+b(k)^2),1)
  d(k) = d(k) * c;
end
%*** change basis
for k = 1:length(Delta)
  h = B(k+1) - B(k);
  c0(k) = nu(k);
  c1(k) = d(k);
  c2(k) =  3/h^2*(nu(k+1)-h*d(k)-nu(k))-1/h  *(d(k+1)-d(k));
  c3(k) = -2/h^3*(nu(k+1)-h*d(k)-nu(k))+1/h^2*(d(k+1)-d(k));
end
c0'
c1'
c2'
c3'
  
figure(1)
plot([0,B(1)],[nu(1),nu(1)],'rs-');
hold on
for k = 1 : length(nu)-1
  s = linspace(B(k),B(k+1),20);
  smz0 = s-B(k);
  plot(s, c0(k)+(c1(k)+(c2(k)+c3(k)*smz0).*smz0).*smz0,'r-');
  plot([B(k),B(k+1)],[nu(k),nu(k+1)],'ks');
end

figure(2)
plot([0,B(1)],[d(1),d(1)],'rs-');
hold on
for k = 1 : length(nu)-1
  s = linspace(B(k),B(k+1),20);
  smz0 = s-B(k);
  plot(s, c1(k)+(2*c2(k)+3*c3(k)*smz0).*smz0,'r-');
  plot([B(k),B(k+1)],[d(k),d(k+1)],'gs');
end

figure(1)
s = linspace(B(end),2.5*B(end),200);
k1 = B(end)   * ( 2 * (nu_inf-nu(end)) - B(end) * d(end) );
k2 = B(end)^2 * ( nu_inf - nu(end) - B(end) * d(end) );
fs = nu_inf + (-k1 + k2./s)./s ;
k1,k2,nu_inf
B'
plot(s,fs,'r-');


plot(B,nu,'k-')
hold off

function [x, N] = cg( A, b, x0)
%% conjugate gradient method
%% author: Hangfeng Gong

%% N iteration times
N = 0;
norm_b = norm(b);

x = x0;
r = b - A * x; %%r0
p = r; %%p0

while 1
    
    if norm(r) / norm_b < 10^(-6)
      break
    end
  
  alpha = r.' * r /(p.' * A * p);
  
  x = x + alpha * p;
  
  N = N + 1;
  
  rr = r - alpha * A * p; 
  
  beta = rr.' * rr / (r.' * r);
  
  p = rr + beta * p;
  
  r = rr;
  
end
   
end


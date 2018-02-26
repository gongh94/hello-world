function [x, N]= sd( A, b, x0 )
%% steepest descent method
%% Copyright (c) 2018 Hangfeng Gong

%% N iteration times

x = x0;
N = 0;

norm_b = norm(b);
 
r = b - A*x;

y_err = norm(r) / norm_b; %% for plotting convergence history
x_N = 0; %% for plotting convergence history
 
    while 1
    
    if norm(r) / norm_b < 10^(-6)
        break
    end
    
    alpha = r.' * r/(r.' * A * r);
    x = x + alpha * r;
    
    N = N+1;
    
    r = b - A*x;
    end
   
end


function [x, N] = pcg( A, b, x0 )
%% Preconditioned Conjugate Gradient Method
%% Copyright (c) 2018 Hangfeng Gong

%% N iteration times

    N = 0;
    
    norm_b = norm(b);
    
    L = ichol(sparse(A)); %% incomplete Cholesky factorization
    
    L = full(L);
    
    M = L * L.';
    
    Minv = inv(M);
    
    x = x0;
    
    r = b - A * x0; %% r0
    
    p = Minv * r; %% p0
    
    z = p; %% z0

    while 1
        
        if norm(r) / norm_b < 10^(-6) 
            break
        end
        
        alpha = r.' * z / (p.' * A * p);
        
        x = x + alpha * p;
        
        N = N + 1;
        
        rr = r - alpha * A * p;
        
        zz = Minv * rr;
        
        beta = rr.' * zz / ( r.' * z);
        
        p = zz + beta * p;
        
        r = rr;
        
        z = zz;
    end
   
end


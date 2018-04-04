function u = exact_soln( t, x)

% x: a row vector; t: a scalar
% u: a col vector

% exact solution to the head equation v_t = niu * v_xx
% v(0,x) = sin(x/2pi)
% v(t,0)= 0
% v(t,1)= 0
    
    niu = 1;
    
    u = exp(- 4* pi *pi *niu * t) * sin( 2*pi*x.');


end
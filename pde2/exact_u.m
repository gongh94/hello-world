function u = exact_u( x, t )
% x: a row vector; t: a scalar
% exact solution to the head equation u_t - u_xx =0

    u = exp(- pi*pi*t) * sin(pi* transpose(x));


end


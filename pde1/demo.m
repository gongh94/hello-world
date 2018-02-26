format long

%% Consider the one-dimension Poisson equation,
%% -d^2(u)/dx^2 = sin(x) over (0,pi) subject to
%% the Dirichlet boundary condition u(0)=0, u(pi)=0.

%% Generate the matrix A and the right-hand side vector b
%% by discretizing the domain and consider second-order finite
%% difference method to construct the discrete Laplace operator Delta_h.

%% To solve Ax = b with A symmetric and positive definite,
%% we will use three methods: the steepest descent method,
%% the conjugate gradient method, and the preconditioned conjugate gradient method.

    n = 10; % n divisions, user-defined
    
    h = pi/n; % length of steps
    
    % construct the matrix A
    v= -1*ones(1,n);
    A = 2*eye(n+1) + diag(v,1)+diag(v,-1); % Ax =b
    
    % note that A has to be symmetric and positive definite 
    A(1,1)=1;
    A(1,2)=0;
    A(2,1)=0; % u(0)=0
    
    A(n+1,n+1)=1;
    A(n+1,n)=0;
    A(n,n+1)=0; % u(pi)=0
    
    % construct the vector b
    b = ones(n+1,1);
    
    for i = 0:n
        b(i+1) = h*h*sin(i*h);
    end
    
    % try some initial guess
    x = ones(n+1,1); 
    %x = (1:1:n+1).'; 
    
    
    [result1, N1] = sd(A,b,x);
    result1
    
    [result2, N2] = cg(A,b,x);
    result2
    
    [result3, N3] = pcg(A,b,x);
    result3

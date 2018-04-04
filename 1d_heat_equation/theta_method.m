function u = theta_method (theta)

% 2018(c) Hangfeng Gong

h = 0.01; % length of spatial step
h_square = h*h;
k = 0.01; % length of temporal step

a = 0;
b = 1;
sigma = 1;
L = b - a;
m = L/h; % m evenly spaced points over space domain

x = linspace(a,b,m);
t = 0;

u = sin(pi * x).'; %initial u
%u(1) = 0;
%u(m) = 0;

aux1 = - sigma / h_square * (1 - theta) * ones(1, m-1);
M1 = (1/k + 2 * sigma / h_square * (1 - theta)) * eye(m) + diag(aux1,1) + diag(aux1,-1);
M1(1,1) = 1;
M1(1,2) = 0;
M1(m,m) = 1;
M1(m,m-1) = 0;

M1_inv = inv(M1);

aux2 = - sigma / h_square * theta * ones(1, m-1);
M2 = ( 2 * sigma / h_square * theta - 1/k) * eye(m) + diag(aux2,1) + diag(aux2,-1);
M2(1,1) = 1;
M2(1,2) = 0;
M2(m,m) = 1;
M2(m,m-1) = 0;

N = 200; 
store = zeros(m, N+1);
store_ = zeros(m, N+1);

store(:,1) = u;
store_(:,1) = exact_u(x,t);
    
for j = 1:N
    t = k * j; % stop at time k*N
    ave_f = (1-theta)*hf(x, t+k) + theta*hf(x,t); % weighted average of f value
    u = M1_inv * (ave_f - M2 * u);
    
    store(:,j+1) = u;
    store_(:,j+1) = exact_u(x,t);
    
end

%u_ = exact_u(x,t);

surf(store);
%surf(store_);

end

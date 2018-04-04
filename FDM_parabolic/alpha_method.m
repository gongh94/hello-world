function [u, eu]= alpha_method (alpha)

% u: col vector; a solution obtained from FDM
% eu: col vector; exact solution

% 2018(c) Hangfeng Gong

M = 30; % M divisions on space; give M + 1 points

dx = 1/M;   % domain of time is [0,1]
            % length of spatial step
          
N = 10; % steps on time

T = 0.5; % length of time

dt = T/N; % length of temporal step

niu = 1;

r = niu * dt/(dx*dx);

x = linspace(0, 1, M+1); 
x = x(2:M); % 1 * (M-1)

u = sin(2*pi*x.'); % u(0,x) excluding spatial boundary points; (M-1) * 1

aux1 = - alpha* r * ones(1, M-2);
Q1 = (1 + 2 * alpha * r)* eye(M-1) + diag(aux1,1) + diag(aux1,-1);

aux2 = (1 - alpha)* r * ones(1, M-2);
Q2 = (1- 2 * r + 2*alpha*r) * eye(M-1) + diag(aux2,1) + diag(aux2,-1);

store = zeros(M-1,N+1);
%estore = zeros(M-1,N+1);

store(:,1) = u;
%estore(:,1) = exact_soln(0,x);

for j = 1:N
    u = Q1\(Q2*u);
    store(:,j+1) = u;
    %estore(:,j+1) = exact_soln(j*dt,x);
end

eu = exact_soln(T, x);

plot = surf(store);
set(plot,'LineStyle','none');
title('T=0.5');
zlabel('Temperature');
xlabel('t');
ylabel('x');

end
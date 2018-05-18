sigma = 0.3;
r = 0.05;
B = 100;
K = 105;
T = 1;

M = 20; %% t steps
N = 40; %% x steps

dt = T/M;
h = B/N;

capM = zeros(N+1);

for i= 2:(N+1)
    capM(i,i-1)=h/6;
    capM(i,i)=2*h/3;
end

capM(1,1)=h/3;

for i=1:N
    capM(i,i+1) = h/6;
end
    
A = zeros(N+1);

for i= 2:(N+1)
    j=i-1;
    A(i,i-1)=-(j*h)*(j*h)*sigma*sigma/2/h + r*(j*h)/2;
    A(i,i)=(j*h)*(j*h)*sigma*sigma/h + r*h;
end

A(1,1) = r/2*h;

for i=1:N
    j=i-1;
    A(i,i+1)=-(j*h)*(j*h)*sigma*sigma/2/h - r*(j*h)/2;
end

A = dt/2*A;

u = zeros(N+1,1);

S = linspace(0,B,N+1).';

w = zeros(N+1,1);

%% initial w = u - g in V_0
w_prev = K * ones(N+1,1) - S; 
w_prev(1) = 0;
w_prev(N+1) = 0;

g = zeros(N+1,1);

% initial g
g_prev = zeros(N+1,1); 
g_prev(1) = K;
g_prev(N+1) = 0; %% singular point on the boundary

% initial u
u_prev = w_prev + g_prev;

% u^{m+1} - u^{m}
Du = zeros(N+1,1);

for m = 1:M
    
    g(1)= K*exp(-r*m*dt);
    g(N+1)= (K-B)*exp(-r*m*dt);

    rhs = -A*(g + g_prev) - capM*(g - g_prev) + (capM-A)*w_prev;
    w = (capM+A)\rhs;
    
    u = w + g;
    
    
    Du = u - u_prev;
    eta2_m = 0;
    
    for k = 1:N
        eta2_m = eta2_m + (Du(k+1)-Du(k))*(Du(k+1)-Du(k))*h/3;
    end
    
    eta2_m = dt * sigma * sigma/2 * eta2_m
    %{
    eta2_mw = 0;
    
   for k = 1:N
       alpha = Du(k)/dt + r * u(k);
       beta = Du(k+1)/dt + r* u(k+1);
       C = r/h*(u(k)-u(k+1));
       fun = @(S) (alpha*k + (beta - alpha)/h*S - beta*(k-1) + C*S)*2;
       q = integral(fun,(k-1)*h,k*h);
       
       eta2_mw = eta2_mw + 1/k*q;
   end
   
   eta2_mw = dt/(sigma*sigma)*eta2_mw
   %}
    
    u_prev = u;
    w_prev = w;
    g_prev = g;
end

u
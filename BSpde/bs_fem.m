% copyright(C)2018 Hangfeng Gong %
%{
price a down-and-in European put via solving the Black-Scholes PDE by using
finite element method;
see more details in the documentation.
%}

S = 90; % stock price; 

sigma = 0.3; % volatility
r = 0.05; % interest rate
B = 100; % barrier
K = 105; % strike 
T = 1; %  expiry

%% uniform time size and space size
nsteps_t = 20; %% #steps on time domain
dt = T/nsteps_t;

nsteps_s = 40; %% #steos on space domain
h = B/nsteps_s;

M = zeros( nsteps_s+1); 

for i= 2:(nsteps_s+1)
    M(i,i-1) = h/6;
    M(i,i) = 2*h/3;
end

M(1,1) = h/3;

for i=1:nsteps_s
    M(i,i+1) = h/6;
end
    
A = zeros(nsteps_s+1);

for i= 2:(nsteps_s+1)
    j=i-1;
    A(i,i-1)=-(j*h)*(j*h)*sigma*sigma/2/h + r*(j*h)/2;
    A(i,i)=(j*h)*(j*h)*sigma*sigma/h + r*h;
end

A(1,1) = r/2*h;

for i=1:nsteps_s
    j=i-1;
    A(i,i+1)=-(j*h)*(j*h)*sigma*sigma/2/h - r*(j*h)/2;
end

A = dt/2*A;

%% w: w(t,x) = u(t,x) - G(t,x) 
%% w_prev: w(t-dt,x)
w = zeros(nsteps_s+1,1);

%% discretize the space domain
x = linspace(0, B, nsteps_s+1).'; 
w_prev = K * ones(nsteps_s+1,1) - x; 
w_prev(1) = 0;
w_prev(nsteps_s+1) = 0;

%% g: G(t,x)
%% g_prev: G(t-dt,x)

g = zeros(nsteps_s+1,1);

g_prev = zeros(nsteps_s+1,1); % initial ug
g_prev(1) = K;
g_prev(nsteps_s+1) = 0; %% singular point

% interpolation for option value at stock price S
H = 2*h;
S_idx = ceil(S/h); %% fine the interval where S falls in

v = ones(1,nsteps_t+1); % u(t,S)
dv = ones(1,nsteps_t+1);

v(1) = w_prev(S_idx);
dv(1) = (w_prev(S_idx+1) - 2*w_prev(S_idx) + w_prev(S_idx-1))/H;

for m = 1:nsteps_t
    
    g(1)= K*exp(-r*m*dt);
    g(nsteps_s+1)= (K-B)*exp(-r*m*dt);

    rhs = -A*(g + g_prev) - M*(g - g_prev) + (M-A)*w_prev;
    w = (M+A)\rhs;
    
    u = w + g; %% u(t,x) 
    
    % store option value and the delta at stock price S
    v(m+1) = u(S_idx);
    dv(m+1) = (u(S_idx+1) - 2*u(S_idx) + u(S_idx-1))/H;
    
    w_prev = w;
    g_prev = g;
end


t = linspace(0,T,nsteps_t+1);
subplot(2,1,1);
plot(t,v)
title('v-t')
xlabel('t')
ylabel('v')

subplot(2,1,2);
plot(t,dv)
title('delta-t')
xlabel('t')
ylabel('delta')

%copyright(C)2018 Hangfeng Gong 
%{
price a down-and-in European put via solving the Black-Scholes PDE by using
finite difference method;
see more details in the documentation.
%}

S = 90;

sigma = 0.3; % volatility
r = 0.05; % interest rate
B = 100; % barrier
K = 105; % strike 
T = 1; %  expiry

nsteps_t = 20;
dt = T/nsteps_t;

nsteps_s = 40;
h = B/nsteps_s;

u0 = zeros(1,nsteps_t+1); %% u(t,0)
uB = zeros(1,nsteps_t+1); %% u(t,B)

for i=1:nsteps_t+1
    %%j=i-1;
    u0(i)=K*exp(-r*(i-1)*dt);
    uB(i)=(K-B)*exp(-r*(i-1)*dt);
end

uB(1)=0; %% singular point

x = linspace(0,B,nsteps_s + 1).';
u = K - x; %% u(0,x)

% interpolation for option value at stock price S
H = 2*h;
S_idx = ceil(S/h); %% fine the interval where S falls in

v = ones(1,nsteps_t+1); % u(t,S)
dv = ones(1,nsteps_t+1);

v(1) = u(S_idx);
dv(1) = (u(S_idx+1) - 2*u(S_idx) + u(S_idx-1))/H;

u = u(2:N); %% u_1,...,u_{N-1}, except u_0 and u_N for each time node
newS_idx = S_idx - 1;

for k = 2:(nsteps_t+1)
    gamma = 1/2*r*dt;
    
    Ml = zeros(N-1,N-1);
    Mr = zeros(N-1,N-1);
    
    for i = 1:(nsteps_s-1)
        alpha = sigma*sigma/4*i*i*dt;
        beta = r*i/4*dt;
        Ml(i,i) = 1+2*alpha+gamma;
        Mr(i,i) = 1-2*alpha-gamma;
        
        if i ~= 1
        Ml(i,i-1)= -alpha + beta;
         Mr(i,i-1) = alpha - beta;
        end
        
        if i ~= nsteps_s-1
            Ml(i,i+1)= -alpha-beta;
            Mr(i,i+1) = alpha+beta;
        end
    end
    
    aux = Mr*u;
    
    alpha1 = sigma*sigma/4*dt;
    beta1 = r/4*dt;
    aux(1) = aux(1) + (alpha1-beta1)*(u0(k)+u0(k-1));
    
    alphaB = sigma*sigma*(N-1)*(N-1)/4*dt;
    betaB = r/4*(N-1)*dt;
    aux(N-1) = aux(N-1) + (alphaB-betaB)*(uB(k)+uB(k-1));
    
    u = Ml\aux;
    
    v(k) = u(newS_idx);
    dv(k) = (u(newS_idx+1) - 2*u(newS_idx) + u(newS_idx-1))/H;
end

%%u = [u0(nsteps_t+1);u;uB(nsteps_t+1)];

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



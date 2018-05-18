format long;

N = 64; %% N intervals, N+1 points

x = linspace(0,1,N+1);

h = 1/N;

aux= -1/h*ones(1,N);
A = 2/h*eye(N+1) + diag(aux,1) + diag(aux,-1);
A(1,1)=1;
A(N+1,N+1)=1;
A(1,2)=0;
A(N+1,N)=0;

c = sqrt(2)/2; 
F = zeros(N+1,1);

for i = 1:N
    if(x(i)<c && x(i+1)<c) %% indicate the interval with singular point of f
        F(i)= h;
    end
    
    if(x(i)<c && x(i+1)>c)
        F(i)= 3/4*h;
        F(i+1)= 1/4*h;
    end
    
    %% other elements of F is zero
end

F(1)=0;

cffns = A\F;

%{
somex = 0.3;

for i=1:N
    if(x(i)< somex && x(i+1) > somex)
        xidx=i;
    end
end

if(xidx == 1)
    ksi = 2*(somex-x(1))/h - 1;
    sim_u = cffns(1)*(1-ksi)/2;
elseif(xidx == N)
    ksi = 2*(somex - x(N))/h -1;
    sim_u = cffns(N+1)*(1+ksi)/2;
else
    ksi = 2*(somex - x(xidx))/h - 1;
    sim_u = cffns(i)*(1-ksi)/2 + cffns(i+1)*(1+ksi)/2;
end

sim_u
exact_u = exact(somex,c)
%}

%cffns(1)

E = 0;

for i= 1:N+1
    e = abs(cffns(i) - exact(x(i),c));
    if(e>E)
        E=e;
    end
end
    
E
    
    
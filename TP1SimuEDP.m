close all
clear all

N=14; % dim espace

M=2*N*N; % dim temps


h = 1/(N-1);
k = 1/(M-1);
alpha = 1;
r=0.5;

xi=[0:h:1];


% CONSTRUCTION DE C

C = -2*eye(N);

for i=1:N

    for j=1:N

        if i == j+1 || i == j-1

            C(i,j) = 1;

        end
    end
end

% CONSTRUCTION DE A

Aexp = eye(N) + (k/(h*h))*C;
Aimp = eye(N);
Ar = r*Aimp + (1-r)*Aexp;


% CONSTRUCTION DE B

Bexp = eye(N);
Bimp = eye(N) - (k/(h*h))*C;
Br = r*Bimp + (1-r)*Bexp;

% CONSTRUCTION DE U

U = zeros(N,M);

% CONDITION INITIALE



for i =1:N
    x = (i-1)*h;
    U(i,1)= (x/(2*alpha))*(x-1);
    
end

figure(13)
plot(xi,U(:,1));
xlabel("espace");
ylabel("u(x,t)");
title("Condition initiale");




% CONSTRUCTION SOLUTION SCHEMA EXPLICITE

Uexp = U;


for j = 1:M-1
    
    Uexp(:,j+1) = Aexp*Uexp(:,j);

    for i = 1:M
   
        Uexp(1,i)=0;
        Uexp(N,i)=0;           % CONDITIONS LIMITES

    end
    
    figure(1)
    hold on
    plot(xi,Uexp(:,j));
end
figure(1)
plot(xi,Uexp(:,M));

xlabel("espace");
ylabel("u(x,t)");
title("Schéma explicite");

% Solution Analytique

M=2*N*N; % dim temps

N=14; %dim espace

h = 1/(N-1);
k = 1/(M-1);

n = 2000;


xi=[0:h:1];

t=[0:k:1];

nk=[1:1:n];

w = nk*pi;

Kd = 2./(w.^3).*(((-1).^nk)-1);

D = diag(Kd);

Uanaly = sin(xi'*w)*D*exp(-(w.^2)'*t);


for i = 1: M

    plot(xi,Uanaly(:,i));

end

figure(30)

plot(Uanaly)
xlabel("espace");
ylabel("u(x,t)");
title("Solution théorique");


% % CONSTRUCTION SOLUTION SCHEMA IMPLICITE

Uimp = U;
invB = inv(Bimp);

for j = 1:M-1

    Uimp(:,j+1) = invB*Uimp(:,j);

    for i = 1:M

        Uimp(1,i)=0;
        Uimp(N,i)=0;           % CONDITIONS LIMITES

    end

    figure(2)
    hold on
    plot(xi,Uimp(:,j));
end
figure(2)
plot(xi,Uimp(:,M));

xlabel("espace");
ylabel("u(x,t)");
title("Schéma implicite");

% CONSTRUCTION SOLUTION SCHEMA C-N

Ucn = U;
invB = inv(Br);

for j = 1:M-1

    Ucn(:,j+1) = invB*Ar*Ucn(:,j);

    for i = 1:M

        Ucn(1,i)=0;
        Ucn(N,i)=0;           % CONDITIONS LIMITES

    end

    figure(3)
    hold on
    plot(xi,Ucn(:,j));
end
figure(3)
plot(xi,Ucn(:,M));

xlabel("espace");
ylabel("u(x,t)");
title("Schéma de Crank-Nicolson");
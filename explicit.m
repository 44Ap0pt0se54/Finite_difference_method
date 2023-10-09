function [Ucn] = explicit(M,N)

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

% CONSTRUCTION SOLUTION SCHEMA EXPLICITE


Ucn = U;
invB = inv(Br);

for j = 1:M-1

    Ucn(:,j+1) = invB*Ar*Ucn(:,j);

    for i = 1:M

        Ucn(1,i)=0;
        Ucn(N,i)=0;           % CONDITIONS LIMITES

    end

end
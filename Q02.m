% Solving u.t = u.xx in [-10,10]
% u(x=10,t) = 0 Dirichlet
% u.x(x=-10, t) = 0 Noyman

clear

% Parameters
L = 20;  % Length of the rod 
dx = 0.1;  % Delta x
Nx = L/dx;  % Steps of x
dt = 0.002;  % Delta t
Nt = 1/dt;  % Steps of t

x = zeros(1, Nx);
t = zeros(1, Nt);
u = zeros(Nx, Nt);

% Initial conditions

for n = 1:Nx+1
    x(n) = (n-1) * dx - L/2;  % Initiate x
    u(n,1) = u0(x(n));  % Initiate u(x,t=0)
end

% Boundry conditions

for j = 1:Nt+1
    u(Nx+1,j) = 0;  % u(x=10,t) = 0 
    t(j) = (j-1)*dt;  % Initiate t
end

% Applay center difference
% dtu(n,j) = (u(n+1,j) - 2*u(n,j) + u(n-1,j))/(dx^2);

A = zeros(Nx+1,Nx+1);
for i = 1:Nx+1
    if i == 1  % because of Noyman
        A(i,i) = -1;
        A(i,i+1) = 1;
    elseif i < Nx+1
        A(i,i-1) = 1;
        A(i,i) = -2;
        A(i,i+1) = 1;
    % i = Nx+1 A(i,:) = 0 because of Dirichlet
    end
end
A = A/(dx^2);

% Solve u't = A*u

% Trapzoidal rule 

for i=2:Nx  % x axis
    for j=2:Nt  % t axis
        d1 = A(i,:)*u(:,j);
        d2= A(i,:)*u(:,j-1);
        u(i,j) = u(i,j-1)+dt*0.5*(d1+d2);
    end
end

% plot the solution in different times

figure('position',[550,250,900,700]);
subplot(2,3,1);
plot(x,u(:,1))
title('t=0')
subplot(2,3,2);
plot(x,u(:,20))
title('t=20')
subplot(2,3,3);
plot(x,u(:,100))
title('t=100')
subplot(2,3,4);
plot(x,u(:,250))
title('t=250')
subplot(2,3,5);
plot(x,u(:,400))
title('t=400')
subplot(2,3,6);
plot(x,u(:,500))
title('t=500')

function u = u0(x)  % define u at t = 0
    if abs(x) <= 1
        u = 1;
    else
        u = 0;
    end
end
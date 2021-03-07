% reformulate the optimal time problem as a fixed end-time problem
close all
clear all
clc

[PAR, CON, SC, SCu] = par_robot;
nx = 4;
nu = 2;
nte = 1;
dtau = 0.01;

N = 100;
T = N*dtau;
tau = 0:dtau:T;

% initial conditions
x0 = [-5+2*pi; 0; -4+2*pi; 0];

% weights
%Q = 100*diag([1,1,1,1]);
R = 1000*diag([1]);

% constraints
umin    =   [CON.u1(1); CON.u2(1)];
umax    =   [CON.u1(2); CON.u2(2)];
xmin    =   [CON.x1(1); CON.x2(1); CON.x3(1); CON.x4(1)];
xmax    =   [CON.x1(2); CON.x2(2); CON.x3(2); CON.x4(2)];

xSS = [pi/2; 0; 0; 0];
uSS = [0; 0];

import casadi.*
ocp = casadi.Opti();

X       =   ocp.variable(nx,N+1);
U       =   ocp.variable(nu,N);
TE      =   ocp.variable(nte,N);
X0      =   ocp.parameter(nx,1);

J       =   0;

for i=1:N
    % dynamics
    xx   =   rk4_freetime(@(t,x,u)robot_ode(t,x,u),dtau,tau(i),X(:,i),U(:,i),TE(:,i));
    ocp.subject_to( X(:,i+1) == xx);
    % constraints
    ocp.subject_to( SCu.*umin <= U(:,i) <= SCu.*umax );
    ocp.subject_to( xmin <= X(:,i+1) <= xmax );
    ocp.subject_to( TE(:,i) > 0 );
    % cost
    dx      =   X(:,i+1) - xSS;
    dte = TE(:,i);
    %J       =   J + 0.5*dx'*Q*dx + 0.5*dte'*R*dte;
    J       =   J + dte'*R*dte;
end

% initial condition constraint
ocp.subject_to( X(:,1) == X0 );
% terminal cost
%J       =   J + 1000*(X(:,end) - xSS)'*(X(:,end) - xSS);
% terminal constraint
ocp.subject_to( X(:,end) == xSS );

% solve initial OCP
% set initial condition
ocp.set_value(X0, x0);

% initial guess
ocp.set_initial(X,repmat(x0,1,N+1))
ocp.set_initial(U,repmat(SCu.*(umax-umin)/2,1,N));
ocp.set_initial(TE,repmat(5,1,N));

% set objective
ocp.minimize(J);
% solve ocp
ocp.solver('ipopt');
sol = ocp.solve();

% get solution
Xsol    =   sol.value(X);
Usol    =   sol.value(U);
TEsol   =   sol.value(TE);

% plot results
figure(1)
subplot(3,1,1)
plot(tau,Xsol, 'Linewidth', 2)
title('state trajectories')

subplot(3,1,2)
plot(tau,[Usol,Usol(:,end)], 'Linewidth', 2)
title('input trajectories')

subplot(3,1,3)
plot(tau,[TEsol,TEsol(:,end)], 'Linewidth', 2)
title('optimal end time')
% closed loop MPC control
close all
clear all
clc


nx = 4;
nu = 2;
dt = 0.03;

Nmpc = 100;
Tmpc = Nmpc*dt;
tmpc = 0:dt:Tmpc;

% simulation time
t0      =   0; %start time
tf      =   3; %end time of simulation
tsim    =   0:dt:tf;
ksim    =   length(tsim);

% initial conditions
x0 = [-5+2*pi; 0; -4+2*pi; 0];

% weights
Q = 100*diag([1,1,1,1]);
R = diag([1,1]);

[PAR, CON, SC, SCu] = par_robot;
% constraints
umin    =   [CON.u1(1); CON.u2(1)];
umax    =   [CON.u1(2); CON.u2(2)];
xmin    =   [CON.x1(1); CON.x2(1); CON.x3(1); CON.x4(1)];
xmax    =   [CON.x1(2); CON.x2(2); CON.x3(2); CON.x4(2)];

xSS = [pi/2; 0; 0; 0];
uSS = [0; 0];

import casadi.*
ocp = casadi.Opti();

X       =   ocp.variable(nx,Nmpc+1);
U       =   ocp.variable(nu,Nmpc);
X0      =   ocp.parameter(nx,1);

J       =   0;
for i=1:Nmpc
    % dynamics
    xx      =   rk4(@(t,x,u)robot_ode(t,x,u),dt,tmpc(i),X(:,i),U(:,i));
    ocp.subject_to( X(:,i+1) == xx);
    % constraints
    ocp.subject_to( SCu.*umin <= U(:,i) <= SCu.*umax );
    ocp.subject_to( xmin <= X(:,i+1) <= xmax );
    % cost
    dx      =   X(:,i+1) - xSS;
    du      =   U(:,i) - uSS;
    J       =   J + 0.5*dx'*Q*dx + 0.5*du'*R*du;
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
ocp.set_initial(X,repmat(x0,1,Nmpc+1))
ocp.set_initial(U,repmat(SCu.*(umax-umin)/2,1,Nmpc));

% set objective
ocp.minimize(J);
% solve ocp
ocp.solver('ipopt');
sol = ocp.solve();

% get solution
Xsol    =   sol.value(X);
Usol    =   sol.value(U);

plot_robot((0:size(Usol,2))*dt,Xsol,Usol,xSS,uSS,[5,6])

%% MPC simulation with wrong parameters
[xsim,usim] = MPCloop_robot(U,X,X0,ocp,@(t,x,u)wrong_robot_ode(t,x,u),dt,tf,x0,SCu.*(umax-umin)/2);
wrong_plot_robot(tsim,xsim,usim,xSS,uSS,[7,8])

%% MPC simulation with Gaussian measurement error
%[xsim,usim] = MPCloop_robot(U,X,X0,ocp,@(t,x,u)robot_ode(t,x,u),dt,tf,x0,SCu.*(umax-umin)/2);
%wrong_plot_robot(tsim,xsim,usim,xSS,uSS,[7,8])

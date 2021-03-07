% open loop optimal control
close all
clear all
clc

nx = 4;
nu = 2;
dt = 0.03;

Nmpc = 100;
Tmpc = Nmpc*dt;
tmpc = 0:dt:Tmpc;

% initial conditions
x0 = [-5+2*pi; 0; -4+2*pi; 0];

% weights
Q = 100*diag([1,1,1,1]);
R = diag([1,1]);

% open loop control with correct parameters
[PAR, CON, SC, SCu] = par_robot;
% constraints
umin    =   [CON.u1(1); CON.u2(1)];
umax    =   [CON.u1(2); CON.u2(2)];
xmin    =   [CON.x1(1); CON.x2(1); CON.x3(1); CON.x4(1)];
xmax    =   [CON.x1(2); CON.x2(2); CON.x3(2); CON.x4(2)];

% set point
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
    
    %xx      =   eulerf(@(t,x,u)robot_ode(t,x,u),dt,tmpc(i),X(:,i),U(:,i));
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
J       =   J + 1000*(X(:,end) - xSS)'*(X(:,end) - xSS);
% terminal constraint
%ocp.subject_to( X(:,end) == xSS );

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

% verify solution via accurate integration
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
t_ode = [];
z_ode = [];

z0i = x0;
for i=1:Nmpc
    t_span = [(i-1)*dt, i*dt];
    [t_i, z_i] = ode45(@robot_ode, t_span, z0i, options, Usol(:,i));
    t_ode = [t_ode; t_i];
    z_ode = [z_ode; z_i];
    z0i = z_ode(end, :);
end

% plot results
plot_robot((0:size(Usol,2))*dt,Xsol,Usol,xSS,uSS,[1,2])
figure(1)
hold on
subplot(4,1,1)
plot(t_ode, z_ode(:,1), 'g', 'Linewidth', 2)
legend('steady state', 'x_1', 'x_1 from ode45', 'Location', 'Southeast')
subplot(4,1,2)
plot(t_ode, z_ode(:,2), 'g', 'Linewidth', 2)
legend('steady state', 'x_2', 'x_2 from ode45', 'Location', 'Southeast')
subplot(4,1,3)
plot(t_ode, z_ode(:,3), 'g', 'Linewidth', 2)
legend('steady state', 'x_3', 'x_3 from ode45', 'Location', 'Southeast')
subplot(4,1,4)
plot(t_ode, z_ode(:,4), 'g', 'Linewidth', 2)
legend('steady state', 'x_4', 'x_4 from ode45', 'Location', 'Southeast')

% plot the results in the x-y plane and q1-q2 plane
l1 = 0.5;
l2 = 0.5;
x1_plot = zeros(1,Nmpc+1);
y1_plot = zeros(1,Nmpc+1);
x2_plot = zeros(1,Nmpc+1);
y2_plot = zeros(1,Nmpc+1);

for i = 1:Nmpc+1
    x1_plot(i) = l1 * cos(Xsol(1,i));
    y1_plot(i) = l1 * sin(Xsol(1,i));
    x2_plot(i) = x1_plot(i) + l2 * cos(Xsol(1,i)+Xsol(3,i));
    y2_plot(i) = y1_plot(i) + l2 * sin(Xsol(1,i)+Xsol(3,i));
end

figure(3)      % plot in the x-y plane
subplot(2,1,1)
plot(x1_plot, y1_plot)
xlabel('x_1')
ylabel('y_1')
title('plot results in the x-y plane')
subplot(2,1,2)
plot(x2_plot, y2_plot)
xlabel('x_2')
ylabel('y_2')

figure(4)      % plot in the q1-q2 plane
plot(Xsol(1,:), Xsol(3,:))
xlabel('q_1')
ylabel('q_2')
title('plot results in the q1-q2 plane')



%% open loop simulation with wrong parameters
[PAR, CON, SC, SCu] = wrong_par_robot;
% constraints
umin    =   [CON.u1(1); CON.u2(1)];
umax    =   [CON.u1(2); CON.u2(2)];
xmin    =   [CON.x1(1); CON.x2(1); CON.x3(1); CON.x4(1)];
xmax    =   [CON.x1(2); CON.x2(2); CON.x3(2); CON.x4(2)];

xSS = [pi/2; 0; 0; 0];
uSS = [0; 0];

x_sim = zeros(nx,Nmpc+1);
u_sim = Usol;
x0j = x0;
x_sim(:,1) = x0; 
for j=1:Nmpc
    % dynamics
    xx_sim      =   rk4(@(t,x,u)wrong_robot_ode(t,x,u),dt,tmpc(j),x0j,u_sim(:,j));
    x_sim(:,j+1) = xx_sim;
    x0j = xx_sim;
    
end

wrong_plot_robot((0:size(Usol,2))*dt,x_sim,u_sim,xSS,uSS,[5,6])

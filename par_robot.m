function [PAR, CON, SC, SCu] = par_robot
% definition of system parameters
PAR.b1 = 200.0;
PAR.b2 = 50.0;
PAR.b3 = 23.5;
PAR.b4 = 25.0;
PAR.b5 = 122.5;
PAR.c1 = -25.0;
PAR.g1 = 784.8;
PAR.g2 = 245.3;
PAR.l1 = 0.5;
PAR.l2 = 0.5;

% constraints
CON.x1  = [-pi, pi]; 
CON.x2  = [(-3/2)*pi,  (3/2)*pi];
CON.x3  = [-pi, pi];
CON.x4  = [(-3/2)*pi,  (3/2)*pi];
CON.u1  = [-1000,  1000];
CON.u2  = [-1000,  1000];

% scaling factors
SC.u2 = 1E-2;
SC.u1 = 1E-2;

SCu     =   [SC.u1; SC.u2];
end

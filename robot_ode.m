function dx = robot_ode(t,x,u)
    [PAR, CON, SC, SCu] = par_robot;
    b1 = PAR.b1;
    b2 = PAR.b2;
    b3 = PAR.b3;
    b4 = PAR.b4;
    b5 = PAR.b5;
    c1 = PAR.c1;
    g1 = PAR.g1*SC.u1;
    g2 = PAR.g2*SC.u1;
    l1 = PAR.l1;
    l2 = PAR.l2;
    
    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    x4 = x(4);
    u1 = u(1);
    u2 = u(2);
    
    dx1 = x2;
    dx3 = x4;
    
    B1 = b1 + b2*cos(x3);
    B2 = b3 + b4*cos(x3);
    B3 = b3 + b4*cos(x3);
    B4 = b5;
    B = [B1, B2; B3, B4];
    B_inv = (1/(B1*B4-B2*B3)).*[B4, -B2; -B3, B1];
    
    dx24 = B_inv * ( c1*sin(x3)*[x2^2+x2*x4+x4^2; -x2^2] -...
                            [g1*cos(x1)/SC.u1+g2*cos(x1+x3)/SC.u1; g2*cos(x1+x3)/SC.u1]...
                            + [u1/SC.u1; u2/SC.u2] );
    dx2 = dx24(1);
    dx4 = dx24(2);
    dx = [dx1; dx2; dx3; dx4];
end


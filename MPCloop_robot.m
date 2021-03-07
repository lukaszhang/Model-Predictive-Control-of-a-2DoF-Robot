function [x,u] = MPCloop_robot(U,X,X0,ocp,myfun,dt,tf,xinit,uinit)
% inputs:
%   X0      optipar for initial condition for each OCP
%   U       optivar for MPC inputs
%   sol     parameterized solution to OCP
%   x0      numerical value of initial condition
    t           =   0:dt:tf;
    N           =   length(t)-1;
    nx          =   length(xinit);
    [nu,Nmpc]   =   size(U);
    x           =   zeros(nx,N+1);
    u           =   zeros(nu,N);
    x(:,1)      =   xinit;
    Xsol        =   repmat(xinit,1,Nmpc+1);
    Usol        =   repmat(uinit,1,Nmpc);
    
    for k = 1:N
        % set current state
        %x(:,k) = x(:,k) + 0.1 .* randn(nx,1);      % add Gaussian measurement error
        ocp.set_value(X0, x(:,k)); 
        % initialize OCP
        ocp.set_initial(X,Xsol)
        ocp.set_initial(U,Usol);
        % resolve OCP
        sol = ocp.solve();
        Xsol    =   sol.value(X);
        Usol    =   sol.value(U);
        % get first element
        u(:,k)   =   Usol(:,1);
        % apply to system using rk4
        x(:,k+1) =   rk4(@(t,x,u)myfun(t,x,u),dt,t(k),x(:,k),u(:,k));

    end
end
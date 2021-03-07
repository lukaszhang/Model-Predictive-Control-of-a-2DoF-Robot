% rk4 integrator for the reformulated free end-time problem
function xf = rk4_freetime(ode,h,t,x,u,te)
  k1 = te.*ode(t,x,u);
  k2 = te.*ode(t,x+h/2*k1,u);
  k3 = te.*ode(t,x+h/2*k2,u);
  k4 = te.*ode(t,x+h*k3,  u);
  xf = x + h/6 * (k1 + 2*k2 + 2*k3 + k4); 
end
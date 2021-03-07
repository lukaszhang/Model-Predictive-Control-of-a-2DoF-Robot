function plot_robot(t,x,u,xSS,uSS,fignum)
linewidth = 1.5;
[~,~,~,SCu] = par_robot; 
ksim    =   length(t);
% rescale results
x    =   x./repmat([1;1;1;1],1,ksim);
u    =   u./repmat(SCu,1,ksim-1);
xSSp    =   xSS;
uSSp    =   uSS./SCu;
warning('The results have to be normalized for plotting.')
figure(fignum(1))
subplot(4,1,1)
plot(t,xSSp(1)*ones(size(t)),'LineWidth',linewidth);
hold on;
plot(t,x(1,:),'Linewidth',linewidth);
grid on;
xlabel('$t$ [s]','interpreter','Latex');
ylabel('$q1(t)$ [rad]','interpreter','Latex');

subplot(4,1,2)
plot(t,xSSp(2)*ones(size(t)),'LineWidth',linewidth);
hold on;
plot(t,x(2,:),'Linewidth',linewidth);
grid on;
xlabel('$t$ [s]','interpreter','Latex');
ylabel('$\dot{q1}(t)$ [rad/s]','interpreter','Latex');

subplot(4,1,3)
plot(t,xSSp(3)*ones(size(t)),'LineWidth',linewidth);
hold on;
plot(t,x(3,:),'Linewidth',linewidth);
grid on;
xlabel('$t$ [s]','interpreter','Latex');
ylabel('$q2(t)$ [rad]','interpreter','Latex');

subplot(4,1,4)
plot(t,xSSp(4)*ones(size(t)),'LineWidth',linewidth);
hold on;
plot(t,x(4,:),'Linewidth',linewidth);
grid on;
xlabel('$t$ [s]','interpreter','Latex');
ylabel('$\dot{q2}(t)$ [rad/s]','interpreter','Latex');

figure(fignum(2))
subplot(2,1,1)
plot(t,uSSp(1)*ones(size(t)),'LineWidth',linewidth);
hold on;
stairs(t(1:end-1),u(1,:),'Linewidth',linewidth);
grid on
xlabel('$t$ [s]','interpreter','Latex');
ylabel('$u_1^\star(t)$ [Nm]','interpreter','Latex');

subplot(2,1,2)
plot(t,uSSp(2)*ones(size(t)),'LineWidth',linewidth);
hold on;
stairs(t(1:end-1),u(2,:),'Linewidth',linewidth);
grid on
xlabel('$t$ [s]','interpreter','Latex');
ylabel('$u_2^\star(t)$ [Nm]','interpreter','Latex');
%%% Check the funtions to obtain Euler Angles
n = length(time);

subplot(3,1,1)
plot(  time(1:n), rad2deg(Euler_ang_true(1,1:n)))
hold on
plot(  time(1:n), rad2deg(Euler_ang_true2(1,1:n)))
ax=gca;
ax.YGrid = 'on';
ax.XLim = [0   time(n)];
xlabel('$ t $ : Time $(s)$','Interpreter', 'latex ')
ylabel('$ \phi$ : Roll Angle $(^{\circ})$','Interpreter','latex')
legend('True by qtoEuler' ,'True by CtoEuler','Interpreter','latex');

subplot(3,1,2)
plot(  time(1:n), rad2deg(Euler_ang_true(2,1:n)))

hold on
plot(  time(1:n), rad2deg(Euler_ang_true2(2,1:n)))
ax=gca;
ax.YGrid = 'on';
ax.XLim = [0   time(n)];
xlabel('$ t $ : Time $(s)$','Interpreter', 'latex ')
ylabel('$ \theta $ : Pitch Angle $(^{\circ})$','Interpreter','latex')
legend('True by qtoEuler' ,'True by CtoEuler','Interpreter','latex');

subplot(3,1,3)
plot(  time(1:n),rad2deg(Euler_ang_true(3,1:n)))
hold on
plot(  time(1:n),rad2deg(Euler_ang_true2(3,1:n)))
ax=gca;
ax.YGrid = 'on';
ax.XLim = [0  time(n)];
xlabel('$ t $ : Time $(s)$','Interpreter', 'latex ')
ylabel('$ \psi$ : Yaw Angle $(^{\circ})$','Interpreter','latex')
legend('True by qtoEuler' ,'True by CtoEuler','Interpreter','latex');

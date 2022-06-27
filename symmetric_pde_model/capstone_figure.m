%% Generate Capstone Figure
simple_advection_1D_exp;

figure(2);

subplot(3,2,1)
plot(x(2:end),rho(:,1),'b'); % Empirical
hold on;
plot(x12,pow2(a*(1*T/nt)).*(log(16))*pow2(-x12),'Color','Red'); % Analytic
hold off;
axis([x_bounds(1) x_bounds(2) 0 6]);
legend('Empirical','Analytic');
title('Exponential IC: t = 0')

subplot(3,2,3)
plot(x(2:end),rho(:,nt/2),'b'); % Empirical
hold on;
plot(x12,pow2(a*((nt/2)*T/nt)).*(log(16))*pow2(-x12),'Color','Red'); % Analytic
hold off;
axis([x_bounds(1) x_bounds(2) 0 6]);
title('Exponential IC: t = 1')

subplot(3,2,5)
plot(x(2:end),rho(:,nt),'b'); % Empirical
hold on;
plot(x12,pow2(a*(nt*T/nt)).*(log(16))*pow2(-x12),'Color','Red'); % Analytic
hold off;
axis([x_bounds(1) x_bounds(2) 0 6]);
title('Exponential IC: t = 2')

clear all;

simple_advection_1D_uni;

figure(2);

subplot(3,2,2)
plot(x(2:end),rho(:,1),'b'); % Empirical
hold on;
plot(x12,pow2(a*(1*T/nt)).*(log(16))*pow2(-x12),'Color','Red'); % Analytic
hold off;
axis([x_bounds(1) x_bounds(2) 0 6]);
legend('Empirical','Analytic');
title('Uniform IC: t = 0')

subplot(3,2,4)
plot(x(2:end),rho(:,nt/2),'b'); % Empirical
hold on;
plot(x12,pow2(a*((nt/2)*T/nt)).*(log(16))*pow2(-x12),'Color','Red'); % Analytic
hold off;
axis([x_bounds(1) x_bounds(2) 0 6]);
title('Uniform IC: t = 1')

subplot(3,2,6)
plot(x(2:end),rho(:,nt),'b'); % Empirical
hold on;
plot(x12,pow2(a*(nt*T/nt)).*(log(16))*pow2(-x12),'Color','Red'); % Analytic
hold off;
axis([x_bounds(1) x_bounds(2) 0 6]);
title('Uniform IC: t = 2')

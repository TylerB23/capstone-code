%% Specify Parameters

lambda = 1; % advection speed (growth rate)
D = 0.0001; % diffusion constant
k = 0.01; % time step size
h = lambda*k; % nutrient step size 
x_bounds = [0 2];
T = 2;

nx = floor((x_bounds(2) - x_bounds(1))/h); % number of nutrient steps
nt = T/k;

x = x_bounds(1):h:x_bounds(2); % grid points
t = 0:k:T;

rho = zeros(nx,nt); % density function

%% Initial Condition

xtmp = x(1<= x & x < 2);
rho(1 <= x & x < 2,1) = pow2(-xtmp/lambda);
rho(1 <= x & x < 2,1) = rho(1 <= x & x < 2,1)/(h*sum(rho(1 <= x & x < 2,1)));

%% Euler

i = 2:nx-1;% free node indices
c1 = (lambda*k/(2*h));
c2 = (lambda*lambda*k*k/(2*h*h));
c3 = sqrt(2*D);

for dt = 1:nt
    
    rho(i,dt+1) = rho(i,dt) - k*lambda*rho(i-1,dt) + c3*randn(198,1);
    % Doubling Boundary Conditions
    rho(x==1,dt+1) = 2*rho(end,dt);
    rho(end,dt+1) = rho(end,dt) - k*lambda*rho(end-1,dt) + c3*randn;
    % Zero below x=1
    rho(x < 1, dt+1) = 0;
end

%% Plotting

rho_max = max(rho, [], 'all');
rho_min = min(rho, [], 'all');
x12 = x(x>=1); % only part of x-axis where analytic sol'n is defined

figure(1);
M(nt) = struct('cdata',[],'colormap',[]); % gathers plot frames for animation

for dt = 1:nt
    plot(x(2:end),rho(:,dt),'b'); % Empirical
    hold on;
    plot(x12,pow2(lambda*(dt*T/nt)).*(log(16))*pow2(-x12),'Color','Red'); % Analytic
    hold off;
    axis([x_bounds(1) x_bounds(2) rho_min 1.1*rho_max]);
    legend('Empirical','Analytic');
    drawnow;
    M(dt) = getframe;
end

% To display animation:
% movie(M,1,10)

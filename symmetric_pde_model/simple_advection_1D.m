%% Specify Parameters

a = 1; % advection speed (growth rate)
h = 0.05; % nutrient step size
k = 0.05; % time step size
x_bounds = [0 2];
T = 2;

nx = (x_bounds(2) - x_bounds(1))/h; % number of nutrient steps
nt = T/k;

x = x_bounds(1):h:x_bounds(2); % grid points
t = 0:k:T;

rho = zeros(nx,nt); % density function

%% Initial Conditions

% uniform distribution around 1
rho(1-5*h <= x & x <= 1+5*h,1) = round(rand(sum(1-5*h <=x & x<= 1+5*h),1));
% cut off stuff below 1
rho(x<=1,1) = 0;
% ensure rho(x==1,1) has a nonzero value
rho(x==1,1) = 1;

%% Lax-Friedrichs

% build matrix
% A = diag(zeros(nx)) + diag(ones(nx),1) + diag(-1*ones(nx),-1);
% A(end:1) = 2;
% A(1:end) = -2;
% A = (-a/(2*h))*A;
% i = 2:nx-1;% free node indices

% for dt=1:nt
%     
%     rho(i,dt+1) = 0.5*(rho(i-1,dt)+rho(i+1,dt)) - (a*k/(2*h))*(rho(i+1,dt)-rho(i-1,dt));
%     % Periodic Boundary Conditions
%     rho(x==1,dt+1) = 0.5*(rho(x==1+h,dt)) - (a*k/(2*h))*(rho(x==1+h,dt) - 2*rho(end,dt));
%     rho(end,dt+1) = 0.5*(rho(end-1,dt)) - (a*k/(2*h))*(rho(x==1,dt)-rho(end-1,dt));
%     
% end

%% Lax-Wendroff

i = 2:nx-1;% free node indices

for dt = 1:nt
    
    rho(i,dt+1) = rho(i,dt) - (a*k/(2*h))*(rho(i+1,dt)-rho(i-1,dt)) + (a*a*k*k/(2*h*h))*(rho(i-1,dt) - 2*rho(i,dt) + rho(i+1,dt));
    % Doubling Boundary Conditions
    rho(x==1,dt+1) = 2*rho(end,dt);
    rho(end,dt+1) = rho(end,dt) - (a*k/(2*h))*(rho(x==1,dt)-rho(end-1,dt)) + (a*a*k*k/(2*h*h))*(rho(end-1,dt) - 2*rho(end,dt) + rho(x==1,dt));
    
end

%% Beam-Warming

% define i to account for one-sided method
% i = 3:nx-1;
% 
% for dt = 1:nt
%     
%     rho(i,dt+1) = rho(i,dt) - (a*k/(2*h))*(3*rho(i,dt)-4*rho(i-1,dt)+rho(i-2,dt)) + ...
%                   (a*a*k*k/(2*h*h))*(rho(i,dt)-2*rho(i-1,dt)+rho(i-2,dt));
%     rho(x==1+h,dt+1) = rho(x==1+h,dt) - (a*k/h)*(rho(x==1+h,dt)-rho(x==1,dt));
%     % Boundary conditions
%     rho(x==1,dt+1) = 2*rho(end,dt);
%     rho(end,dt+1) = rho(end,dt) - (a*k/(2*h))*(3*rho(end,dt)-4*rho(end-1,dt)+rho(end-2,dt)) + ...
%                     (a*a*k*k/(2*h*h))*(rho(end,dt)-2*rho(end-1,dt)+rho(end-2,dt));
% 
% end

%% Plotting

rho_max = max(rho, [], 'all');
rho_min = min(rho, [], 'all');

figure(1);
M(nt) = struct('cdata',[],'colormap',[]); % gathers plot frames for animation

for dt = 1:nt
    plot(x(2:end),rho(:,dt),'b'); % Empirical
    hold on;
    plot(x,pow2(a*(dt*T/nt)).*rho(x==1,1)*pow2(-x),'Color','Red'); % Analytic
    hold off;
    axis([x_bounds(1) x_bounds(2) rho_min 1.1*rho_max]);
    legend('Empirical','Analytic');
    drawnow;
    M(dt) = getframe;
end

% To display animation:
% movie(M,1,10)

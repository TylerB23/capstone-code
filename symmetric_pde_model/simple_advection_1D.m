%% Specify Parameters

a = 1; % advection speed
h = 0.005; % spatial step size
k = 0.005; % time step size
x_bounds = [0 2];
time_bounds = [0 2];

nx = (x_bounds(2) - x_bounds(1))/h; % number of x steps
nt = (time_bounds(2) - time_bounds(1))/k;

x = x_bounds(1):h:x_bounds(2); % grid points

rho = zeros(nx,nt); % density function

%% Initial Conditions

% rho(x==1,1) = 1;

% uniform distribution around 1
rho(1-5*h <= x & x <= 1+5*h,1) = round(rand(sum(1-5*h <=x & x<= 1+5*h),1));
% cut off stuff below 1
rho(x<=1,1) = 0;

%% Boundary Conditions

% free node indices
i = 2:nx-1;

%% Lax-Friedrichs

% build matrix
% A = diag(zeros(nx)) + diag(ones(nx),1) + diag(-1*ones(nx),-1);
% A(end:1) = 2;
% A(1:end) = -2;
% A = (-a/(2*h))*A;

% for dt=1:nt
%     
%     rho(i,dt+1) = 0.5*(rho(i-1,dt)+rho(i+1,dt)) - (a*k/(2*h))*(rho(i+1,dt)-rho(i-1,dt));
%     % Periodic Boundary Conditions
%     rho(x==1,dt+1) = 0.5*(rho(x==1+h,dt)) - (a*k/(2*h))*(rho(x==1+h,dt) - 2*rho(end,dt));
%     rho(end,dt+1) = 0.5*(rho(end-1,dt)) - (a*k/(2*h))*(rho(x==1,dt)-rho(end-1,dt));
%     
% end

%% Lax-Wendroff

% for dt = 1:nt
%     
%     rho(i,dt+1) = rho(i,dt) - (a*k/(2*h))*(rho(i+1,dt)-rho(i-1,dt)) + (a*a*k*k/(2*h*h))*(rho(i-1,dt) - 2*rho(i,dt) + rho(i+1,dt));
%     % Doubling Boundary Conditions
%     % rho(x==1,dt+1) = rho(x==1,dt) - (a*k/(2*h))*(rho(x==1+h,dt)-rho(end,dt)) + (a*a*k*k/(2*h*h))*(rho(end,dt) - 2*rho(x==1,dt) + rho(x==1+h,dt));
%     rho(x==1,dt+1) = 2*rho(end,dt);
%     rho(end,dt+1) = rho(end,dt) - (a*k/(2*h))*(rho(x==1,dt)-rho(end-1,dt)) + (a*a*k*k/(2*h*h))*(rho(end-1,dt) - 2*rho(end,dt) + rho(x==1,dt));
%     
% end

%% Beam-Warming

% define i to account for one-sided method
i = 3:nx-1;

for dt = 1:nt
    
    rho(i,dt+1) = rho(i,dt) - (a*k/(2*h))*(3*rho(i,dt)-4*rho(i-1,dt)+rho(i-2,dt)) + ...
                  (a*a*k*k/(2*h*h))*(rho(i,dt)-2*rho(i-1,dt)+rho(i-2,dt));
    rho(x==1+h,dt+1) = rho(x==1+h,dt) - (a*k/h)*(rho(x==1+h,dt)-rho(x==1,dt));
% Boundary conditions
%   rho(x==1,dt+1) = rho(x==1,dt) - (a*k/(2*h))*(3*rho(x==1,dt)-2*4*rho(end,dt)+rho(end-1,dt)) + ...
%                    (a*a*k*k/(2*h*h))*(rho(x==1,dt)-2*rho(end,dt)+rho(end-1,dt));
    rho(x==1,dt+1) = 2*rho(end,dt);
    rho(end,dt+1) = rho(end,dt) - (a*k/(2*h))*(3*rho(end,dt)-4*rho(end-1,dt)+rho(end-2,dt)) + ...
                  (a*a*k*k/(2*h*h))*(rho(end,dt)-2*rho(end-1,dt)+rho(end-2,dt));

end

%% Plotting

rho_max = max(rho, [], 'all');
rho_min = min(rho, [], 'all');
figure(1);

% put time in top of graph

for dt = 1:nt
    plot(x(2:end),rho(:,dt),'b');
    axis([x_bounds(1) x_bounds(2) rho_min rho_max]);
    %axis([0 2 0 2]);
    pause
end

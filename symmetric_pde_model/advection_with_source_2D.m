%%
% advection without source term in two spatial dimensions

%% Specifying Parameters

a = [1 2]; % velocity vector - may want to make into a function later
h = 0.1; % spatial step size
k = 0.05; % time step size
x_bounds = [0 3];
y_bounds = [0 3];
t_bounds = [0 1];
nx = (x_bounds(2) - x_bounds(1))/h; % number of x steps
ny = (y_bounds(2) - y_bounds(1))/h; % number of y steps
nt = (t_bounds(2) - t_bounds(1))/k; % number of time steps

x = x_bounds(1):h:x_bounds(2); % grid points
y = y_bounds(1):h:y_bounds(2); % grid points

rho = zeros(nx,ny,nt); % initialize density matrix

%% Initial Conditions

% simplest case
rho(x==1,y==1,1) = 1;

% uniform distribution in area 1 < x,y < 2
% the `sum` calls add up to the number of entries which satisfy the
% conditions so the LHS and RHS have the same dimensions
%
% rho(1<=x & x<=2,1<=y & y<=2,1) = round(rand(sum(1<=x & x<=2),sum(1<=y & y<=2)));

%% Boundary conditions
% should set these after solving at each point in time, since I need to do
% that anyway as the values may get overwritten

rho(x==0,y==0,:) = 0;

% free node indices
i = 2:nx-1;
j = 2:ny-1;

%% Lax-Wendroff Method

for dt=1:nt
    
    rho(i,j,dt+1) = rho(i,j,dt) - ...
                    (a(1)*k/(2*h))*(rho(i+1,j,dt)-rho(i-1,j,dt)) - ...
                    (a(2)*k/(2*h))*(rho(i,j+1,dt)-rho(i,j-1,dt)) + ...
                    (a(1)*a(1)*k*k/(h*h*2))*(rho(i+1,j,dt)-2*rho(i,j,dt)+rho(i-1,j,dt)) + ...
                    (a(2)*a(2)*k*k/(h*h*2))*(rho(i,j+1,dt)-2*rho(i,j,dt)+rho(i,j-1,dt));
    % Boundary Conditions
    % rho(x==1,y>=1,dt+1) = 
    % rho(x>=1,y==1,dt+1) = 
    % rho(x==2,y>=2,dt+1) = ...;
    % rho(x>=2,y==2,dt+1) = ...;
    % or...
    rho(x>=1,y>=1,dt+1) = rho(x>=1,y>=1,dt) - ...
                          (a(1)*k/(2*h))*(rho(x>=1+h,y>=1,dt)-rho(x>=2,y>=1,dt)) - ...
                          (a(2)*k/(2*h))*(rho(x>=1,y>=1+h,dt)-rho(x>=1,y>=2,dt)) + ...
                          (a(1)*a(1)*k*k/(h*h*2))*(rho(x>=1+h,y>=1,dt)-2*rho(x>=1,y>=1,dt)+rho(x>=2,y>=1,dt)) + ...
                          (a(2)*a(2)*k*k/(h*h*2))*(rho(x>=1,y>=1+h,dt)-2*rho(x>=1,y>=1,dt)+rho(x>=1,y>=2,dt));
    rho(x>=2,y>=2,dt+1) = rho(x>=2,y>=2,dt) - ...
                          (a(1)*k/(2*h))*(rho(x>=1,y>=2,dt)-rho(x>=2-h,y>=2,dt)) - ...
                          (a(2)*k/(2*h))*(rho(x>=2,y>=1,dt)-rho(x>=2,y>=2-h,dt)) + ...
                          (a(1)*a(1)*k*k/(h*h*2))*(rho(x>=1,y>=2,dt)-2*rho(x>=2,y>=2,dt)+rho(x>=2-h,y>=2,dt)) + ...
                          (a(2)*a(2)*k*k/(h*h*2))*(rho(x>=2,y>=1,dt)-2*rho(x>=2,y>=2,dt)+rho(x>=2,y>=2-h,dt));
end

% Upwind method
% 
% for dt=1:nt
%     rho(i,j,dt+1) = rho(i,j,dt) - ...
%                     (a(1)*k/h)*(rho(i,j,dt)-rho(i-1,j,dt)) - ...
%                     (a(2)*k/h)*(rho(i,j,dt)-rho(i,j-1,dt)) + ... 
%                     Second order terms below aren't changed from LW yet
%                     (a(1)*a(1)*k*k/(h*h*2))*(rho(i+1,j,dt)-2*rho(i,j,dt)+rho(i-1,j,dt)) + ...
%                     (a(2)*a(2)*k*k/(h*h*2))*(rho(i,j+1,dt)-2*rho(i,j,dt)+rho(i,j-1,dt));
% end


%% Plotting

for dt=1:nt
    hHM = heatmap(rho(:,:,dt))
    hHM.NodeChildren(3).YDir='normal'; % turn Y-Axis normal direction
    pause
end
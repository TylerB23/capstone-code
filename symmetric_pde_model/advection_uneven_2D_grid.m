%% Advection without source in 2 dimensions on an uneven grid

%% Parameters and Grid

a = [0.5 1]; % particle velocity
h = 0.1;  % spatial step
k = 0.05;  % time step

T = 10;    % total time for simulation

% the three parts of the grid are adjusted to make the middle section 
% finer and to remove overlaps
x = [0:h:1 (1+h/2):(h/2):2 (2+h):h:4];
y = [0:h:1 (1+h/2):(h/2):2 (2+h):h:4];

nx = length(x);
ny = length(y);
nt = T/k;

rho = zeros(nx,ny,nt);

%% Initial Conditions

% simplest case
rho(x==1,y==1,1) = 1;

% uniform distribution in area 1 < x,y < 2
% the `sum` calls add up to the number of entries which satisfy the
% conditions, so the LHS and RHS have the same dimensions

% rho(1<=x & x<=2,1<=y & y<=2,1) = round(rand(sum(1<=x & x<=2),sum(1<=y & y<=2)));

%% Beam-warming method

i = 3:nx-1; % x index variable
j = 3:ny-1; % y index variable

% logical indices for boundary lines are:

for dt=1:nt
    
    % Step One: Advection
    rho(i,j,dt+1) = rho(i,j,dt) - ...
                    (a(1)*k/(2*h))*(3*rho(i,j,dt)-4*rho(i-1,j,dt)+rho(i-2,j,dt)) - ...
                    (a(1)*k/(2*h))*(3*rho(i,j,dt)-4*rho(i,j-1,dt)+rho(i,j-2,dt)) + ...
                    (a(1)*a(1)*k*k/(h*h*2))*(rho(i,j,dt)-2*rho(i-1,j,dt)+rho(i-2,j,dt)) + ...
                    (a(1)*a(1)*k*k/(h*h*2))*(rho(i,j,dt)-2*rho(i,j-1,dt)+rho(i,j-2,dt));
    % Step Two: Handle Splitting Boundaries
    rho(x==1,(2>=y)&(y>=1),dt) = 2*rho(x==2,(4>=y)&(y>=2),dt);
    % rho(x==2,(4>=y)&(y>=2),dt) = 0;
    rho((2>=x)&(x>=1),y==1,dt) = 2*rho((4>=x)&(x>=2),y==2,dt);
    % rho((4>=x)&(x>=2),y==1,dt)= 0;
    rho(4>=x & x>=2, 4>=y & y>=2,dt) = 0;
end

%% Visualize

rho_max = max(rho,[],'all');
rho_min = min(rho,[],'all');

%% Heatmap
% for dt=1:nt
%     H = heatmap(rho(:,:,dt))
%     H.NodeChildren(3).YDir='normal'; % turn Y-Axis normal direction
%     pause
% end

%% surface Visualization
% for dt=1:nt
%     surfc(x,y,rho(:,:,dt));
%     axis([0 4 0 4])
%     title("Time : " + num2str(dt))
%     pause
% end

%% Movie
h = figure(1);
M(T) = struct('cdata',[],'colormap',[]);
h.Visible = 'off';

for dt = 1:nt
    surfc(x,y,rho(:,:,dt));
    axis([0 4 0 4])
    title("Time : " + num2str(dt))    
    
    drawnow
    M(dt) = getframe;
end
h.Visible = 'on';

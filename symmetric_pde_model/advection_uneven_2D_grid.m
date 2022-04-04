%% Advection in 2 dimensions on an uneven grid

%% Parameters and Grid

a = [0.05 0.1]; % particle velocity
% abs_a = sqrt(a(1)*a(1) + a(2)*a(2));
% max_a = max(a);
sum_a = sum(a);

k = 0.1;  % time step
h = 10*sum_a*k; % spatial step
T = 1;    % total time for simulation

% the three parts of the grid are adjusted to make the middle section 
% finer and to remove overlaps
x = [0:h:1-h 1:(h/2):2-h/2 2:h:4+h];
y = [0:h:1-h 1:(h/2):2-h/2 2:h:4+h];

nx = length(x);
ny = length(y);
nt = T/k;

rho = zeros(nx,ny,nt);

%% Initial Conditions

% uniform distribution in area 1 < x,y < 2
% the `sum` calls add up to the number of entries which satisfy the
% conditions, so the LHS and RHS have the same dimensions

rho(1.25<=x & x<=1.75,1.25<=y & y<=1.75,1) = rand(sum(1.25<=x & x<=1.75),sum(1.25<=y & y<=1.75));

%% Beam-warming method

i = 3:nx-1; % x index variable
j = 3:ny-1; % y index variable

for dt=1:nt
    
    % Step One: Advection
    rho(i,j,dt+1) = rho(i,j,dt) - ...
                    (a(1)*k/(2*h))*(3*rho(i,j,dt)-4*rho(i-1,j,dt)+rho(i-2,j,dt)) - ...
                    (a(2)*k/(2*h))*(3*rho(i,j,dt)-4*rho(i,j-1,dt)+rho(i,j-2,dt)) + ...
                    (a(1)*a(1)*k*k/(h*h*2))*(rho(i,j,dt)-2*rho(i-1,j,dt)+rho(i-2,j,dt)) + ...
                    (a(2)*a(2)*k*k/(h*h*2))*(rho(i,j,dt)-2*rho(i,j-1,dt)+rho(i,j-2,dt));
    % Step Two: Handle Splitting Boundaries
    rho(x==1,(2>=y)&(y>=1),dt) = rho(x==2,(4>=y)&(y>=2),dt);
    rho((2>=x)&(x>=1),y==1,dt) = rho((4>=x)&(x>=2),y==2,dt);
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

%% Movie

% swap to heatmap to match paper
% make axes from 0 to 100 in each

g = figure(1);
% M(nt) = struct('cdata',[],'colormap',[]);
% g.Visible = 'off';

for dt = 1:nt
  % surf(x,y,rho(:,:,dt));
    H = heatmap(rho(:,:,dt));
    H.NodeChildren(3).YDir='normal'; % turn Y-Axis normal direction
    % axis([0 4 0 4 rho_min rho_max])
    title("Time : " + num2str(dt))
    
  drawnow;
  % M(dt) = H;
  pause
end
% g.Visible = 'on';

%% frame-by-frame Visualization
% for dt=1:nt
%     surf(x,y,rho(:,:,dt));
%     axis([0 4 0 4 rho_min rho_max])
%     title("Time : " + num2str(dt))
%     pause
% end

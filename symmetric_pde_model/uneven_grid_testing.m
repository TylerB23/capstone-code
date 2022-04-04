%% Make an uneven grid

h = 0.1; % spatial step

% x and y should be twice as dense for values between 1 and 2 as it is for
% values greater than 2

% x12 = linspace(1,2,(2/h)+1);
% y12 = linspace(1,2,(2/h)+1);
% x23 = linspace(2+h,4,2*((1/h)+1));
% y23 = linspace(2+h,4,2*((1/h)+1));

x01 = 0:h:1-h;
y01 = 0:h:1-h;
x12 = 1:h/2:2-h/2;
y12 = 1:h/2:2-h/2;
x23 = 2:h:4+h;
y23 = 2:h:4+h;

x = [x01 x12 x23];
y = [y01 y12 y23];

% visualize
figure;

subplot(2,2,1);
title('x coordinate spacing');
plot(x,x,'o')

subplot(2,2,2);
title('y coordinate spacing');
plot(y,y,'o')

%% Make a matrix with this spacing

nx = length(x);
ny = length(y);

M = zeros(nx,ny);

% check sizes
x_1 = M(x==1, (2>=y)&(y>=1));
x_2 = M(x==2, (4>=y)&(y>=2));
y_1 = M((2>=x)&(x>=1), y==1);
y_2 = M((4>=x)&(x>=2), y==2);

% visualize
M(x==1, (2>=y)&(y>=1)) = 1;
M(x==2, (4>=y)&(y>=2)) = -1;
M((2>=x)&(x>=1), y==1) = 2;
M((4>=x)&(x>=2), y==2) = -2;
subplot(2,2,3);
spy(M)

% surface plot
subplot(2,2,4);
surf(x,y,M)

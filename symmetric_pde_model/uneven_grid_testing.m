%% Make an uneven grid

h = 0.05; % spatial step

% x and y should be twice as dense for values between 1 and 2 as it is for
% values greater than 2

x12 = linspace(1,2,(2/h)+1);
y12 = linspace(1,2,(2/h)+1);

% x23 = linspace(2+h,4,2*((1/h)+1));
% y23 = linspace(2+h,4,2*((1/h)+1));
x23 = (2+h):h:4;
y23 = (2+h):h:4;

x = [x12 x23];
y = [y12 y23];

% visualize
figure;

subplot(2,2,1);
plot(x,x,'o')

subplot(2,2,2);
plot(y,y,'o')

%% Make a matrix with this spacing

nx = length(x);
ny = length(y);

M = zeros(nx,ny);

% check sizes
x_1 = M(x==1, (2>=y)&(y>=1));
x_2 = M(x==2, (4>=y)&(y>=2));

% visualize
M(x==1, (2>=y)&(y>=1)) = 1;
M(x==2, (4>=y)&(y>=2)) = -1;
subplot(2,2,3);
spy(M);

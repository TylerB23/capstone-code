function [particles,counts] = Basic_Sim(dt, T, c, xbar, ybar)

% dt is the time step
% T is the total simulation time
% c is the rate of accumulation of y, dy/dt

% Generate a random initial condition
% each particle is a row like [x y] where each are in (1,2)
n = round(100*rand());
particles = 4*rand(n,2);
counts = [n];

% Prepare figure containing pop. distributions
figure(1)
xlabel("X Nutrient Count")
ylabel("Y Nutrient Count")
title("Population distributions over time")
axis([1, 4 1 4])
hold on

for t = 0:dt:T
    % take the step for each particle
    step = [1*dt,c*dt];
    particles = particles + step;
    
    % check for splitting
    for k = 1:1:length(particles)-1
        if particles(k,1) > xbar && particles(k,2) > ybar
            new_cell = particles(k,:) / 2;
            particles(k,:) = new_cell;
            particles = [particles; new_cell];
        end
    end 
    figure(1)
    counts = [counts;length(particles)];
    plot(particles(:,1),particles(:,2),'.'); % fix graphing - try points, not lines
    hold on
end

% plot counts
figure(2)
xlabel("Time")
ylabel("Count")
title("Particle Count over Time")
plot((1:length(counts))*dt,counts);

end
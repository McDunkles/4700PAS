
% Define start and end points (Say these are in meters)
X1 = 0;
X2 = 10;
res = 0.1;
num_points = ((X2-X1)/res) + 1;


% Wave propagation speed
c = 3 * 10^8;


% Simulation step time, in seconds
dt = 0.5 * (10^-9);

z_coor = 0;


% 1D spatial domain
x = linspace(X1,X2, num_points);

% Starting time is 0
t = 0;

% Starting waves
% Ef = zeros(2, num_points);
% Er = zeros(2, num_points);

Ef = zeros(2, num_points);
for i = 1:10
    Ef(1,i) = i*(10-i);
end


fig = figure();
ax = axes('parent', fig);
plt = plot3(ax, x, Ef(1,:), Ef(2,:));

num_steps = 100;


n = 1.5;

% Steps the simulation
for i = 1:num_steps
    vg = c/n;

    % dz-dt Synchronization, dz given in meters
    dz = vg*dt;

    z_coor = z_coor + dz;

    % step_int = 0;
    % Ef(1,1 + step_int:num_steps+1) = Ef(1,1:num_steps + 1 - step_int);

    Ef(1,2:num_steps+1) = Ef(1,1:num_steps);
    Ef(1,1) = 0;
    set(plt, 'XData', x, 'YData', Ef(1,:), 'ZData',Ef(2,:));
%{
    if (mod(i, 10) == 0)
        i
        Ef
    end
%}

    count = num_steps;
    pause(0.10);
end




cla(ax);
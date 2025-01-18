
% Define start and end points (Say these are in meters)
X1 = 0;
X2 = 10;
res = 0.1;
num_points = ((X2-X1)/res) + 1;


% Wave propagation speed
c = 3 * 10^8;
n = 1.25;
vg = c/n;


% Simulation step time, in seconds
dt = 0.5 * (10^-9);

% dz-dt Synchronization, dz given in meters
dz = vg*dt;
frac_z = 0;


% Reflection Coefficients
RL = [0, 0.9];
RR = [0, 0.9];


% 1D spatial domain
x = linspace(X1,X2, num_points);

% Starting time is 0
t = 0;


Ef = zeros(2, num_points);
Er = zeros(2, num_points);

% Input stream
E_in = zeros(2, num_points);

for i = 1:10
    E_in(1, i) = i*(10-i);
end


fig = figure();
subplot(3, 1, 1);
% hold on 
plt1 = plot(x, Ef(1, :), "-", x, Ef(2, :), "--");

subplot(3, 1, 2);
plt2 = plot(x, Er(1, :), "-", x, Er(2, :), "--");


% Left/Right I/O Graph
subplot(3, 1, 3);
plt3 = plot(x, zeros(1, num_points));

num_steps = 100;

first_fwd = 0;
first_rev = 0;

E_in

total_shift = 0;

% Steps the simulation
for i = 1:10*num_steps

    shift_fl = dz/res;
    shift = cast(shift_fl, "int32");
    frac_z = frac_z + (shift_fl - double(shift));
    
    if (frac_z >= 0.5)
        shift = shift + 1;
        frac_z = frac_z - 1;
    end

    total_shift = total_shift + shift;

    if (i < 20)

        msg_format = "[%d] Shift = %d, frac_z = %f, total_shift = %d\n";
        message = sprintf(msg_format, i, shift, frac_z, total_shift);

        fprintf(message);
    end
    

    % step_int = 0;
    % Ef(1 + step_int:num_steps+1) = Ef(1:num_steps + 1 - step_int);

    last_fwd_entry = Ef(:,num_steps+1);
    last_rev_entry = Er(:,1);

    % Step forward field
    Ef(1, 2:num_steps+1) = Ef(1, 1:num_steps);
    Ef(2, 2:num_steps+1) = Ef(2, 1:num_steps);

    if i <= num_steps
        Ef(1, 1) = E_in(1, i);
        Ef(2, 1) = E_in(2, i);

    else
        Ef(1, 1) = 0;
        Ef(2, 1) = 0;
    end
    

    fwd_real = [RL(1), Er(2,1); RL(2), Er(1,1)];
    fwd_im = [RL(1), -Er(1,1); RL(2), Er(2,1)];

    %{
    if (det(fwd_real) ~= 0) && (first_fwd == 0)
        fwd_real
        det(fwd_real)
        first_fwd = 1
    end
    %}

    Ef(1, 1) = Ef(1, 1) + det(fwd_real);
    Ef(2, 1) = Ef(2, 1) + det(fwd_im);

    set(plt1(1), 'XData', x, 'YData', Ef(1, :));
    set(plt1(2), 'XData', x, 'YData', Ef(2, :));


    % Step reverse field
    Er(1, 1:num_steps) = Er(1, 2:num_steps+1);
    Er(2, 1:num_steps) = Er(2, 2:num_steps+1);

    rev_real = [RR(1), Ef(2, num_steps+1); RR(2), Ef(1, num_steps+1)];
    rev_im = [RR(1), -Ef(1, num_steps+1); RR(2), Ef(2, num_steps+1)];

    Er(1, num_steps+1) = det(rev_real);
    Er(2, num_steps+1) = det(rev_im);

    if (det(rev_im) ~= 0) && (first_rev == 0)
        rev_im
        det(rev_im)
        first_rev = 1
    end

    set(plt2(1), 'XData', x, 'YData', Er(1, :));
    set(plt2(2), 'XData', x, 'YData', Er(2, :));

    count = num_steps;
    pause(0.02);
end


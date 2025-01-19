
% Define start and end points (Say these are in meters)
X1 = 0;
X2 = 10;
res = 0.1;
set_length = ((X2-X1)/res) + 1;


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
x = linspace(X1,X2,set_length);

% Starting time is 0
t = 0;


Ef = zeros(2, set_length);
Er = zeros(2, set_length);

% Input stream
E_in = zeros(2, set_length);

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
plt3 = plot(x, zeros(1, set_length));

num_steps = 100;

first_fwd = 0;
first_rev = 0;

E_in

total_shift = 0;
prev_total_shift = 0;

% Steps the simulation
for i = 1:10*num_steps

    % Compute shift and fraction values
    shift_fl = dz/res;
    shift = cast(shift_fl, "int32");
    frac_z = frac_z + (shift_fl - double(shift));
    
    if (frac_z >= 0.5)
        shift = shift + 1;
        frac_z = frac_z - 1;
    end

    % Update total shift
    total_shift = total_shift + shift;

    % Compliment of shift, from the other side of the array
    % Examples:
    % shift = 1 -> shift_comp = 101
    % shift = 2 -> shift_comp = 100
    % shift_comp = base - shift + 1 = 101 - shift + 1 = 102 - shift
    shift_comp = set_length - shift + 1;

    % Debug message printing shift/frac information
    if (i < 20)

        msg_format = "[%d] Shift = %d, frac_z = %f, total_shift = %d, prev_total_shift = %d\n";
        message = sprintf(msg_format, i, shift, frac_z, total_shift, prev_total_shift);

        fprintf(message);
    end

    
    % ===== Step forward field ===== %

    last_fwd_seg = [Ef(1, shift_comp:set_length); ...
        Ef(2, shift_comp:set_length)];

    % Translate all elements over by 'shift' units
    Ef(1, 1+shift:set_length) = Ef(1, 1:set_length-shift);
    Ef(2, 1+shift:set_length) = Ef(2, 1:set_length-shift);

    % Zero out the outlying first 'shift' rows, ready to be written with
    % new data
    Ef(1, 1:shift) = flip(Er(1, 1:shift), 2);
    Ef(2, 1:shift) = flip(Er(2, 1:shift), 2);
    

    fwd_matrix = [RL(1), -RL(2); RL(2), RL(1)];
    rev_matrix = [RR(1), -RR(2); RR(2), RR(1)];


    fwd_real_seg = reshape(Ef(1, 1:shift), [shift, 1]);
    fwd_im_seg = reshape(Ef(2, 1:shift), [shift, 1]);

    fwd_seg = [fwd_real_seg, fwd_im_seg];

    fwd_new_seg = fwd_seg * fwd_matrix;

    Ef(1, 1:shift) = fwd_new_seg(:,1);
    Ef(2, 1:shift) = fwd_new_seg(:,2);


    if i < 20
        msg_format2 = "E_in[%d] = %d\n";
        message2 = sprintf(msg_format2, i, E_in(1, prev_total_shift+1));

        % E_in(1, prev_total_shift+1:total_shift)

        % Ef(1, 1:total_shift)

        % fprintf(message2);
    end 


    % If the input function is still going, add it to the beginning of Ef
    % Should be filling the first 'shift' elements
    if total_shift <= num_steps
        Ef(1, 1:shift) = Ef(1, 1:shift) + flip(E_in(1, prev_total_shift+1:total_shift), 2);
        Ef(2, 1:shift) = Ef(2, 1:shift) + flip(E_in(2, prev_total_shift+1:total_shift), 2);
    end


    set(plt1(1), 'XData', x, 'YData', Ef(1, :));
    set(plt1(2), 'XData', x, 'YData', Ef(2, :));


    % Step reverse field
    Er(1, 1:set_length-shift) = Er(1, 1+shift:set_length);
    Er(2, 1:set_length-shift) = Er(2, 1+shift:set_length);

    Er(1, shift_comp:set_length) = flip(last_fwd_seg(1, :), 2);
    Er(2, shift_comp:set_length) = flip(last_fwd_seg(2, :), 2);

    % size(last_fwd_seg)
    test = last_fwd_seg(1, :);


    rev_real_seg = reshape(Er(1, shift_comp:set_length), [shift, 1]);
    rev_im_seg = reshape(Er(2,shift_comp:set_length), [shift, 1]);

    rev_seg = [rev_real_seg, rev_im_seg];

    rev_new_seg = rev_seg * rev_matrix;

    Er(1, shift_comp:set_length) = rev_new_seg(:,1);
    Er(2, shift_comp:set_length) = rev_new_seg(:,2);


    %{
    if (det(rev_im) ~= 0) && (first_rev == 0)
        rev_im
        det(rev_im)
        first_rev = 1
    end
    %}


    set(plt2(1), 'XData', x, 'YData', Er(1, :));
    set(plt2(2), 'XData', x, 'YData', Er(2, :));

    % Nonzero values in the imaginary component of the reverse wave
    %{
    rev_nonzero_im = nonzeros(Er(2, :));

    if ((size(rev_nonzero_im, 1) > 0) && (size(rev_nonzero_im, 1) < 9))
        fprintf("==== [%d]==== ", i);
        shift
        rev_nonzero_im

        rev_seg
        rev_matrix
        test
        last_fwd_seg

    end
    %}
    

    pause(0.02);

prev_total_shift = total_shift;
end


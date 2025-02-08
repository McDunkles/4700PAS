% Define start and end points (Say these are in meters)
total_z = 2;
dzi = 0.01;
set_length = ((total_z)/dzi) + 1;





% Number of ticks in the simulation
sim_ticks = 5*(set_length-1);




% Reflection Coefficients
RL = 0.9*(1j);
RR = 0.9*(1j);


% 1D spatial domain
x = linspace(X1,X2,set_length);


% Forward and reverse wave components
Ef = zeros(1, set_length);
Er = zeros(1, set_length);

E_net = Ef + Er;

%
B_env = 0.01 + 0.003*(1j);
% B_env = 0;
B = B_env.*ones(1, set_length);


% Gaussian params
u_gaus = 10;
s_gaus = 4;
amp_gaus = 3;

% Input stream
input_length = 2*u_gaus;
E_in = zeros(1, input_length);
E_in(1, 1:input_length) = f_gaus(amp_gaus, u_gaus, s_gaus, 0:input_length-1)

x_max = 2*max(E_in);

% Now we create the figure that will hold the graphs
fig_graph = figure("Name","Graphs");

% Top Plot - Holds the forward wave component
ax = subplot(3, 2, 1);
plt1 = plot(x, real(Ef(1, :)), "-", x, imag(Ef(1, :)), "--");
title("Forward Component");
ylim([-x_max, x_max]);
xlabel("z (m)");
ylabel("E_f(z)");


% Middle Plot - Holds the forward wave component
subplot(3, 2, 3);
plt2 = plot(x, real(Er(1, :)), "-", x, imag(Er(1, :)), "--");
title("Reverse Component");
ylim([-x_max, x_max]);
xlabel("z (m)");
ylabel("E_r(z)");



% Bottom Plot - Wave Input/Output Amplitudes at the boundaries
subplot(3, 2, 5);

fft_len = 512;

E_net = Ef + Er;
E_fft_in = fft(E_, fft_len);


% This will be the independent variable for the bottom plot
time_dom = (0:sim_ticks-1);

left_input = zeros(1, sim_ticks);
left_output = zeros(1, sim_ticks);
right_input = zeros(1, sim_ticks);
right_output = zeros(1, sim_ticks);

plt3 = plot( ...
    time_dom, left_input, "-r", ...
    time_dom, left_output, "-g", ...
    time_dom, right_input, "-b", ...
    time_dom, right_output, "-m");

title("Waveform at Boundaries");
xlabel("t (ps)");
ylabel("E(z)");

lgd = legend("Left Input", "Left Output", "Right Input", "Right Output");
lgd.FontSize = 6;


% Top right plot - Left Input vs Right Output
subplot(3, 2, 2);

plt4 = plot( ...
    time_dom, left_input, "-r", ...
    time_dom, right_output, "-g");

title("Left Input vs Right Output Waveforms");
xlabel("t (ps)");
ylabel("Right Output");
ylim([-x_max, x_max]);


% Middle right plot - E-Field FFT
subplot(3, 2, 4);
temp_dom = 1:set_length;

plt5 = plot( ...
    temp_dom, left_input, "-r", ...
    temp_dom, right_output, "-b");



title("Left Input vs Right Output Waveforms");
xlabel("t (ps)");
ylabel("Right Output");
ylim([-x_max, x_max]);

lgd = legend("Left Input", "Right Output");
lgd.FontSize = 5;


pause_freq = 3;
break_freq = 10;


% Main simulation
for iter = 1:sim_ticks
    
    % ===== (1) Step forward field ===== %

    % Save the last part of the segment for the reverse field
    last_fwd_seg = Ef(1, set_length);

    % Iterate the forward field component
    % (Ef[t, z] = Ef[t-1,z-1])
    Ef(1, 2:set_length) = Ef(1, 1:set_length-1).*exp(-(1j).*B(2:set_length));
    % Ef(1, 2:set_length) = Ef(1, 1:set_length-1);

    % Update the left side of the forward travelling wave 
    Ef(1, 1) = E_in(1, input_length) + Er(1, 1) * RL;

    set(plt1(1), 'XData', x, 'YData', real(Ef));
    set(plt1(2), 'XData', x, 'YData', imag(Ef));



    % ===== (2) Update graphs of fields at the boundaries =====
    left_input(1, iter) = E_in(1, input_length);
    left_output(1, iter) = abs(RL.*Er(1, 1));
    right_output(1, iter) = abs(RR.*last_fwd_seg);

    set(plt3(1), 'XData', time_dom, 'YData', left_input);
    set(plt3(2), 'XData', time_dom, 'YData', left_output);
    set(plt3(4), 'XData', time_dom, 'YData', right_output);

    % Update top right graph
    set(plt4(1), 'XData', time_dom, 'YData', left_input);
    set(plt4(2), 'XData', time_dom, 'YData', right_output);


    % ===== (3) Step reverse field =====
    Er(1, 1:set_length-1) = Er(1, 2:set_length).*exp(-(1j).*B(2:set_length));
    Er(1, set_length) = last_fwd_seg * RR;

    set(plt2(1), 'XData', x, 'YData', real(Er(1, :)));
    set(plt2(2), 'XData', x, 'YData', imag(Er(1, :)));



    % (4) Finally, step the input function
    E_in(1, 2:input_length) = E_in(1, 1:(input_length-1));
    E_in(1, 1) = 0;
    

    if (mod(iter, break_freq) == 0)
        % k = waitforbuttonpress;
    end


    if (mod(iter, pause_freq) == 0)
        pause(0.001);
    end

end



% Generate gaussian pdf
function func_gaussian = f_gaus(amp_gaus, u_gaus, s_gaus, t)
    % u_gaus = 8;
    % s_gaus = 4;
    func_gaussian = amp_gaus.*(1./(s_gaus.*sqrt(2.*pi))).*exp( -0.5.*((t - u_gaus)./s_gaus).^2 );
end


% Define start and end points (Say these are in meters)
len_z = 0.04;
dzi = 0.0001;
set_length = ((len_z)/dzi) + 1;

dt = 10^-12;

% Number of ticks in the simulation
sim_ticks = 2*(set_length-1);

v_g = dzi/dt;


% Reflection Coefficients
RL = 0.9*(1j);
RR = 0.9*(1j);


% 1D spatial domain
x = linspace(0,len_z,set_length);


% Forward and reverse wave components
Ef = zeros(1, set_length);
Er = zeros(1, set_length);

E_net = Ef + Er;

%
B_env = 0.5 + 0.0001*(1j);
% B_env = 0;
B = B_env.*ones(1, set_length);


% Gaussian params
u_gaus = 100;
s_gaus = 25;
amp_gaus = 6;

% Input stream
input_length = 2*u_gaus;
E_in = zeros(1, input_length);

% fft/input stuff
fft_len = 8000;
sample_rate = 2000;
freq = 500;
period = 1/freq;
omega = 2*pi*freq;
total_time = fft_len/sample_rate;
E_wave = amp_gaus*cos(omega*(0:input_length-1)/sample_rate);

freq_precision = sample_rate/fft_len;

% E_in(1, 1:input_length) = E_wave;
E_in(1, 1:input_length) = f_gaus(amp_gaus, u_gaus, s_gaus, 0:input_length-1)

y_max = 1.5*max(E_in);

% Now we create the figure that will hold the graphs
fig_graph = figure("Name","Graphs");

% Top Plot - Holds the forward wave component
ax = subplot(3, 2, 1);
plt1 = plot(x, real(Ef(1, :)), "-", x, imag(Ef(1, :)), "--");
title("Forward Component");
ylim([-y_max, y_max]);
xlabel("z (m)");
ylabel("E_f(z)");


% Middle Plot - Holds the reverse wave component
subplot(3, 2, 3);
plt3 = plot(x, real(Er(1, :)), "-", x, imag(Er(1, :)), "--");
title("Reverse Component");
ylim([-y_max, y_max]);
xlabel("z (m)");
ylabel("E_r(z)");



% Bottom Plot - Wave Input/Output Amplitudes at the boundaries
subplot(3, 2, 5);


% This will be the independent variable for the bottom plot
time_dom = (0:sim_ticks-1);

left_input = zeros(1, sim_ticks);
left_output = zeros(1, sim_ticks);
right_input = zeros(1, sim_ticks);
right_output = zeros(1, sim_ticks);

plt5 = plot( ...
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

plt2 = plot( ...
    time_dom, left_input, "-r", ...
    time_dom, real(right_output), "-g", ...
    time_dom, imag(right_output), "-b");

title("Left Input vs Right Output Waveforms");
xlabel("t (ps)");
ylabel("Right Output");
ylim([-y_max, y_max]);


% Middle right plot - E-Field FFT

subplot(3, 2, 4);
temp_dom = 1:set_length;

fft_domain = freq_precision*(-fft_len/2:(fft_len/2)-1);

% E_fft_in = fft(E_in, fft_len);
E_fft_in = fft(left_input, fft_len);
fft_in_abs = abs(E_fft_in);
fft_in_abs = fftshift(fft_in_abs);


E_fft_out = fft(right_output, fft_len);
fft_out_abs = abs(E_fft_out);
fft_out_abs = fftshift(fft_out_abs);

plt4 = plot( ...
    fft_domain, fft_in_abs, "-r", ...
    fft_domain, fft_out_abs, "-b");



title("Left Input vs Right Output Waveforms");
xlabel("Frequency (Hz)");
ylabel("Right Output");

fft_min = min(min(fft_in_abs, fft_out_abs));
fft_max = max(max(fft_in_abs, fft_out_abs));

if min(fft_in_abs) == max(fft_in_abs)
    fft_min = fft_min - 1;
    fft_max = fft_max + 1;
end

% ylim([min(fft_in_abs), max(fft_in_abs)]);
ylim([-10, 10]);

% lgd = legend("Left Input", "Right Output");
% lgd.FontSize = 5;



% Bottom right plot - Phase shift

subplot(3, 2, 6);
temp_dom = 1:set_length;

temporal_phase_shift = zeros(1,sim_ticks);

plt6 = plot(time_dom, temporal_phase_shift, "-r");

title("Temporal Phase Shift");
xlabel("Time (ps)");
ylabel("Phase (E_{f})");



pause_freq = 3;
break_freq = 10;

total_phase_shift = 0;
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

    total_phase_shift = total_phase_shift - real(B(iter));
    temporal_phase_shift(1, iter) = total_phase_shift;
    set(plt6(1), 'XData', time_dom, 'YData', temporal_phase_shift);
    


    % ===== (2) Update graphs of fields at the boundaries =====
    left_input(1, iter) = E_in(1, input_length);
    left_output(1, iter) = abs(RL.*Er(1, 1));
    right_output(1, iter) = RR.*last_fwd_seg;

    set(plt5(1), 'XData', time_dom, 'YData', left_input);
    set(plt5(2), 'XData', time_dom, 'YData', left_output);
    set(plt5(4), 'XData', time_dom, 'YData', abs(right_output));

    % Update top right graph
    set(plt2(1), 'XData', time_dom, 'YData', left_input);
    set(plt2(2), 'XData', time_dom, 'YData', real(right_output));
    set(plt2(3), 'XData', time_dom, 'YData', imag(right_output));

    % Update fft plot
    E_fft_in = fft(left_input, fft_len);
    fft_in_abs = abs(E_fft_in);
    fft_in_abs = fftshift(fft_in_abs);
    
    
    E_fft_out = fft(right_output, fft_len);
    fft_out_abs = abs(E_fft_out);
    fft_out_abs = fftshift(fft_out_abs);

    set(plt4(1), 'XData', fft_domain, 'YData', fft_in_abs);
    set(plt4(2), 'XData', fft_domain, 'YData', fft_out_abs);
    


    % ===== (3) Step reverse field =====
    Er(1, 1:set_length-1) = Er(1, 2:set_length).*exp(-(1j).*B(2:set_length));
    Er(1, set_length) = last_fwd_seg * RR;

    set(plt3(1), 'XData', x, 'YData', real(Er(1, :)));
    set(plt3(2), 'XData', x, 'YData', imag(Er(1, :)));



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


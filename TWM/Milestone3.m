% Define start and end points (Say these are in micrometers)
len_z = 1000;
dzi = 1;
set_length = ((len_z)/dzi) + 1;

dz = dzi*10^-6;

% Time ticks are 12.5 femtoseconds
dt = 1.25*10^-14;

% Number of ticks in the simulation
sim_ticks = 2*(set_length-1);



Kappa0 = zeros(1, set_length);
Kappa_start = 1/3;
Kappa_end = 2/3;

kp_s_ind = round(Kappa_start*set_length);
kp_end_ind = round(Kappa_end*set_length);
Kappa0(kp_s_ind:kp_end_ind) = 10000;



% 1D spatial domain
x = linspace(0,len_z,set_length);


% Forward and reverse wave components
Ef = zeros(1, set_length);
Er = zeros(1, set_length);


%
B_env = 0 + 0*(1j);
B = B_env.*ones(1, set_length);

% Multiplies with the marching equations
B_factor = exp(-(1j).*B);


% Gaussian params
u_gaus = 200;
s_gaus = 60;
amp_gaus = 200;

% Input stream
input_length = u_gaus + s_gaus*5;

amp_trig = 20;
freq = 100;
sample_rate = 10000;
% input_trig = f_trig(amp_trig, freq)



E_in_static = zeros(1, input_length);
E_in_static(1, 1:input_length) = f_gaus(amp_gaus, u_gaus, s_gaus, 0:input_length-1);

E_in = E_in_static(1, 1:input_length);

MAX_IN_VAL = max(E_in);


y_max = 1.2*max(E_in);

% Now we create the figure that will hold the graphs
fig_graph = figure("Name","Graphs");

% ===== (1) Top Left Plot - Holds the forward wave component =====
ax = subplot(3, 2, 1);
plt1 = plot(x, real(Ef(1, :)), "-", x, imag(Ef(1, :)), "--");
title("Forward Component");
xlabel("z (m)");
ylabel("E_f(z)");
xticks([200, 400, 600, 800, 1000])
ylim([-y_max/2, y_max]);

% =====


% ===== (3) Middle Left Plot - Holds the reverse wave component =====
subplot(3, 2, 3);
plt3 = plot(x, real(Er(1, :)), "-", x, imag(Er(1, :)), "--");
title("Reverse Component");
xlabel("z (m)");
ylabel("E_r(z)");
ylim([-y_max/2, y_max]);



% ===== (5) Bottom Left Plot - Wave Amplitudes at the boundaries =====
subplot(3, 2, 5);

% This will be the independent variable for the bottom plot
time_domain = dt*(0:sim_ticks-1);

left_input = zeros(1, sim_ticks);
left_output = zeros(1, sim_ticks);
right_input = zeros(1, sim_ticks);
right_output = zeros(1, sim_ticks);

plt5 = plot( ...
    time_domain, left_input, "--m", ...
    time_domain, left_output, "-r", ...
    time_domain, right_input, "--c", ...
    time_domain, right_output, "-b");

title("Waveform at Boundaries");
xlabel("t (ps)");
ylabel("E(z)");

% lgd = legend("Left Input", "Left Output", "Right Input", "Output");
% lgd.FontSize = 6;


% ===== (2) Top Right Plot - Left vs Right Output =====
subplot(3, 2, 2);

plt2 = plot( ...
    time_domain, real(left_output), "-r", ...
    time_domain, imag(left_output), "--m", ...
    time_domain, real(right_output), "-b", ...
    time_domain, imag(right_output), "--c");

title("Output Waves");
xlabel("t (ps)");
ylabel("Output (E)");
ylim([-y_max/2, y_max]);




% ===== (4) Middle Right Plot - E-Field FFT =====

subplot(3, 2, 4);

% fft specifications
fft_len = 5000;
freq_precision = sample_rate/fft_len;
total_time = fft_len/sample_rate;

% Frequencies represented in the fft
fft_domain = freq_precision*(-fft_len/2:(fft_len/2)-1);

% Select how much of the fft to graph, basically selecting the "zoom",
% centered at the zero frequency
view_hw = 0.015;
fft_lb = int32((0.5-view_hw)*fft_len);
fft_ub = int32((0.5+view_hw)*fft_len);

% fft of the input waveform
E_fft_in = fftshift(fft(left_input, fft_len));

% fft of the left and right outputs
E_fft_left = fftshift(fft(left_output, fft_len));
E_fft_right = fftshift(fft(right_output, fft_len));

% Left and right fft - absolute values
fft_left_abs = abs(E_fft_left);
fft_right_abs = abs(E_fft_right);

% These are what actually get plotted
fft_plot1 = fft_left_abs;
fft_plot2 = fft_right_abs;


plt4 = plot( ...
    fft_domain(fft_lb:fft_ub), fft_plot1(fft_lb:fft_ub ), "-r", ...
    fft_domain(fft_lb:fft_ub), fft_plot2(fft_lb:fft_ub ), "-b");

title("Frequency Spectrum");
xlabel("Frequency (Hz)");
ylabel("E");

fft_min = min(min(abs(fft_plot1), abs(fft_plot2)));
fft_max = max(max(abs(fft_plot1), abs(fft_plot2)));

if fft_min == fft_max
    fft_min = fft_min - 1;
    fft_max = fft_max + 1;
end

% ylim([min(fft_in_abs), max(fft_in_abs)]);
ylim([-20, 220]);



% ===== (6) Bottom Right Plot - Phase Plot =====
subplot(3, 2, 6);

% Epsilon value to handle division by zero
epsilon = 10^-8;

% Bode plot associated with the fft graph
% phase_plot1 = atan(imag(fft_plot1)./(real(fft_plot1) + epsilon));
% phase_plot2 = atan(imag(fft_plot2)./(real(fft_plot2) + epsilon ) );
phase_plot1 = angle(E_fft_left);
phase_plot2 = angle(E_fft_right);

% phase_plot1


plt6 = plot( ...
    fft_domain(fft_lb:fft_ub), phase_plot1(fft_lb:fft_ub), "-r", ...
    fft_domain(fft_lb:fft_ub), phase_plot2(fft_lb:fft_ub), "-b");


title("Phase Shift");
xlabel("Frequency (Hz)");
ylabel("Phase (E)");

% ylim([min(fft_in_abs), max(fft_in_abs)]);
ylim([-pi, pi]);


% Supplementary graphs!!
fig_supp = figure("Name","MORE Graphs");

% ===== (7) Top Supp Plot - Input wave fft =====
ax2 = subplot(3, 2, 1);

fft_in_real = real(E_fft_in);
fft_in_imag = imag(E_fft_in);
fft_in_abs = abs(E_fft_in);

plt7 = plot( ...
    fft_domain(fft_lb:fft_ub), fft_in_real(fft_lb:fft_ub), "--r", ...
    fft_domain(fft_lb:fft_ub), fft_in_imag(fft_lb:fft_ub), "--b", ...
    fft_domain(fft_lb:fft_ub), fft_in_abs(fft_lb:fft_ub), "-k");

title("Input FFT - Real and Imag Components");
xlabel("Frequency (Hz)");
ylabel("Input (E_{in})");
ylim([-200, 200]);

% ===== (9) Middle Supp Plot - Holds the reverse wave fft =====
subplot(3, 2, 3);

fft_left_real = real(E_fft_left);
fft_left_imag = imag(E_fft_left);
fft_left_abs = abs(E_fft_left);

plt9 = plot( ...
    fft_domain(fft_lb:fft_ub), fft_left_real(fft_lb:fft_ub), "--r", ...
    fft_domain(fft_lb:fft_ub), fft_left_imag(fft_lb:fft_ub), "--b", ...
    fft_domain(fft_lb:fft_ub), fft_left_abs(fft_lb:fft_ub), "-k");

title("Left FFT - Real and Imag Components");
xlabel("Frequency (Hz)");
ylabel("Left Output (Er)");
ylim([-200, 200]);
% =====

% ===== (11) Bottom Supp Plot - Holds the right output fft =====
subplot(3, 2, 5);

fft_right_real = real(E_fft_right);
fft_right_imag = imag(E_fft_right);
fft_right_abs = abs(E_fft_right);

plt11 = plot( ...
    fft_domain(fft_lb:fft_ub), fft_right_real(fft_lb:fft_ub), "-b", ...
    fft_domain(fft_lb:fft_ub), fft_right_imag(fft_lb:fft_ub), "--c", ...
    fft_domain(fft_lb:fft_ub), fft_right_abs(fft_lb:fft_ub), "-k");

title("Right FFT - Real and Imag Components");
xlabel("Frequency (Hz)");
ylabel("Right Output (Ef)");
ylim([-200, 200]);
% =====

% ===== (10) Middle Right Supp Plot - E_left/E_in TF =====
subplot(3, 2, 4);

% centered at the zero frequency
view_hw_tf = 0.0075;
tf_lb = int32((0.5-view_hw_tf)*fft_len);
tf_ub = int32((0.5+view_hw_tf)*fft_len);


% Left "transfer function"
E_fft_left_tf = E_fft_left./E_fft_in;

fft_left_tf_real = real(E_fft_left_tf);
fft_left_tf_imag = imag(E_fft_left_tf);
fft_left_tf_abs = abs(E_fft_left_tf);

plt10 = plot( ...
    fft_domain(tf_lb:tf_ub), fft_left_tf_real(tf_lb:tf_ub), "--r", ...
    fft_domain(tf_lb:tf_ub), fft_left_tf_imag(tf_lb:tf_ub), "--b", ...
    fft_domain(tf_lb:tf_ub), fft_left_tf_abs(tf_lb:tf_ub), "-k");

title("Left Output TF");
xlabel("Frequency (Hz)");
ylabel("Left Output TF");
ylim([-5, 5]);
% =====


% ===== (12) Bottom Right Supp Plot - E_right/E_in TF =====
subplot(3, 2, 6);

% Right "transfer function"
E_fft_right_tf = E_fft_right./E_fft_in;

fft_right_tf_real = real(E_fft_right_tf);
fft_right_tf_imag = imag(E_fft_right_tf);
fft_right_tf_abs = abs(E_fft_right_tf);

plt12 = plot( ...
    fft_domain(tf_lb:tf_ub), fft_right_tf_real(tf_lb:tf_ub), "--r", ...
    fft_domain(tf_lb:tf_ub), fft_right_tf_imag(tf_lb:tf_ub), "--b", ...
    fft_domain(tf_lb:tf_ub), fft_right_tf_abs(tf_lb:tf_ub), "-k");

title("Right Output TF");
xlabel("Frequency (Hz)");
ylabel("Right Output TF");
ylim([-5, 5]);
% =====


pause_freq = 10;
break_freq = 10;
check_limits_freq = 5;



% Main simulation
for iter = 1:sim_ticks
    
    % ===== (1) Step forward field ===== %

    % Save the last part of the segment for the reverse field
    last_fwd_seg = Ef(1, set_length);

    % Iterate the forward field component
    Ef(1, 2:set_length) = Ef(1, 1:set_length-1).*B_factor(1:set_length-1);

    Ef(1, 1:set_length-1) = Ef(1, 1:set_length-1) + ...
    Er(1, 2:set_length).*(1j).*Kappa0(2:set_length).*dz;

    % Update the left side of the forward travelling wave
    Ef(1, 1) = E_in(1, input_length);


    set(plt1(1), 'XData', x, 'YData', real(Ef));
    set(plt1(2), 'XData', x, 'YData', imag(Ef));




    % ===== (2) Update graphs of fields at the boundaries =====
    left_input(1, iter) = E_in(1, input_length);
    left_output(1, iter) = Er(1, 1);
    right_output(1, iter) = last_fwd_seg;

    left_input_env = real(left_input) + imag(left_input);
    left_output_env = real(left_output) + imag(left_output);
    right_output_env = real(right_output) + imag(right_output);

    set(plt5(1), 'XData', time_domain, 'YData', left_input_env);
    set(plt5(2), 'XData', time_domain, 'YData', left_output_env);
    set(plt5(4), 'XData', time_domain, 'YData', right_output_env);

    % Update top right graph
    set(plt2(1), 'XData', time_domain, 'YData', real(left_output));
    set(plt2(2), 'XData', time_domain, 'YData', imag(left_output));
    set(plt2(3), 'XData', time_domain, 'YData', real(right_output));
    set(plt2(4), 'XData', time_domain, 'YData', imag(right_output));
    


    % Update FFT

    % Update top right graph
    E_fft_in = fftshift(fft(left_input, fft_len));
    E_fft_left = fftshift(fft(left_output, fft_len));
    E_fft_right = fftshift(fft(right_output, fft_len));

    % Absolute values
    fft_left_abs = abs(E_fft_left);
    fft_right_abs = abs(E_fft_right);

    
    fft_plot1 = E_fft_left;
    fft_plot2 = E_fft_right;

    set(plt4(1), 'XData', fft_domain(fft_lb:fft_ub), 'YData', abs(fft_plot1(fft_lb:fft_ub)));
    set(plt4(2), 'XData', fft_domain(fft_lb:fft_ub), 'YData', abs(fft_plot2(fft_lb:fft_ub)));


    % Update phase
    phase_plot1 = angle(E_fft_left);
    phase_plot2 = angle(E_fft_right);

    set(plt6(1), 'XData', fft_domain(fft_lb:fft_ub), 'YData', phase_plot1(fft_lb:fft_ub));
    set(plt6(2), 'XData', fft_domain(fft_lb:fft_ub), 'YData', phase_plot2(fft_lb:fft_ub));

    % Update the supplementary fft graph

    % (7) Top Supplementary Plot: E_in FFT
    fft_in_real = real(E_fft_in);
    fft_in_imag = imag(E_fft_in);
    fft_in_abs = abs(E_fft_in);

    set(plt7(1), 'XData', fft_domain(fft_lb:fft_ub), 'YData', fft_in_real(fft_lb:fft_ub));
    set(plt7(2), 'XData', fft_domain(fft_lb:fft_ub), 'YData', fft_in_imag(fft_lb:fft_ub));
    set(plt7(3), 'XData', fft_domain(fft_lb:fft_ub), 'YData', fft_in_abs(fft_lb:fft_ub));

    % (9) Middle Supplementary Plot: E_left FFT
    fft_left_real = real(E_fft_left);
    fft_left_imag = imag(E_fft_left);
    fft_left_abs = abs(E_fft_left);
    
    set(plt9(1), 'XData', fft_domain(fft_lb:fft_ub), 'YData', fft_left_real(fft_lb:fft_ub));
    set(plt9(2), 'XData', fft_domain(fft_lb:fft_ub), 'YData', fft_left_imag(fft_lb:fft_ub));
    set(plt9(3), 'XData', fft_domain(fft_lb:fft_ub), 'YData', fft_left_abs(fft_lb:fft_ub));


    % (11) Bottom Left Supplementary Plot: E_right FFT
    fft_right_real = real(E_fft_right);
    fft_right_imag = imag(E_fft_right);
    fft_right_abs = abs(E_fft_right);

    set(plt11(1), 'XData', fft_domain(fft_lb:fft_ub), 'YData', fft_right_real(fft_lb:fft_ub));
    set(plt11(2), 'XData', fft_domain(fft_lb:fft_ub), 'YData', fft_right_imag(fft_lb:fft_ub));
    set(plt11(3), 'XData', fft_domain(fft_lb:fft_ub), 'YData', fft_right_abs(fft_lb:fft_ub));


    % (10) Middle Right Supplementary Plot: E_left Transfer Function
    % Now let's take output over input
    E_fft_left_tf = E_fft_left./E_fft_in;

    fft_left_tf_real = real(E_fft_left_tf);
    fft_left_tf_imag = imag(E_fft_left_tf);
    fft_left_tf_abs = abs(E_fft_left_tf);

    set(plt10(1), 'XData', fft_domain(tf_lb:tf_ub), 'YData', fft_left_tf_real(tf_lb:tf_ub));
    set(plt10(2), 'XData', fft_domain(tf_lb:tf_ub), 'YData', fft_left_tf_imag(tf_lb:tf_ub));
    set(plt10(3), 'XData', fft_domain(tf_lb:tf_ub), 'YData', fft_left_tf_abs(tf_lb:tf_ub));


    % (12) Bottom Right Supplementary Plot: E_right Transfer Function

    E_fft_right_tf = E_fft_right./E_fft_in;

    fft_right_tf_real = real(E_fft_right_tf);
    fft_right_tf_imag = imag(E_fft_right_tf);
    fft_right_tf_abs = abs(E_fft_right_tf);

    set(plt12(1), 'XData', fft_domain(tf_lb:tf_ub), 'YData', fft_right_tf_real(tf_lb:tf_ub));
    set(plt12(2), 'XData', fft_domain(tf_lb:tf_ub), 'YData', fft_right_tf_imag(tf_lb:tf_ub));
    set(plt12(3), 'XData', fft_domain(tf_lb:tf_ub), 'YData', fft_right_tf_abs(tf_lb:tf_ub));


    % ===== (3) Step reverse field =====
    Er(1, 1:set_length-1) = Er(1, 2:set_length).*B_factor(2:set_length);
    
    Er(1, 2:set_length) = Er(1, 2:set_length) + ...
    Ef(1, 1:set_length-1).*(1j).*Kappa0(1:set_length-1).*dz;


    set(plt3(1), 'XData', x, 'YData', real(Er(1, :)));
    set(plt3(2), 'XData', x, 'YData', imag(Er(1, :)));



    % (4) Finally, step the input function
    E_in(1, 2:input_length) = E_in(1, 1:(input_length-1));
    E_in(1, 1) = 0;


    if (Ef(1, 201) == MAX_IN_VAL)
        iter
        k = waitforbuttonpress;
    end

    if (mod(iter, check_limits_freq) == 0)
    
    end
    

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

% Generate trig function
function trig_func = f_trig(amp, freq, t)
    trig_func = amp.*e^(-(1j)*2*pi*freq*t);
end



% ========== Domain Specification ==========

% Simulation Domain
ds = 1;
sim_ub = 1000;
Ns = ((sim_ub)/ds) + 1;
sim_domain = linspace(0,sim_ub,Ns);

% Spatial step length (dz/ds)
z_len = 10^-3;
dz = 10^-6;
Nz = 1001;
z_domain = linspace(0, z_len, Nz);

% Time Domain
vp0 = 1.0*10^8;
dt = dz/vp0; % Synchronization condition

t_sim = 2.0*10^-11;
Nt = int32(t_sim/dt) + 1;
time_domain = linspace(0, t_sim, Nt);



% ========== Forward and Reverse Wave+Polarization Components ==========

% Forward and reverse wave components
Ef = zeros(1, Ns); % Ef(n) exists at point zn, time tn
Er = zeros(1, Ns); % Er(n) exists at point zn+0.5dz, time tn+0.5dt

% Forward and reverse polarizations
Pf = zeros(1, Ns);
Pr = zeros(1, Ns);


% Spatial derivative matrices
Def = -speye(Ns)/dz;
Der = speye(Ns)/dz;

for ii=(1:Ns)
    Def(ii,ii+1) = 1/dz;

    if ii > 1
        Der(ii,ii-1) = -1/dz;
    end
end


% Previous fields
Efp = Ef;
Erp = Er;
Pfp = Pf;
Prp = Pr;


% Carrier frequency
xc = 40000;
kc = 2*pi*xc;
% f_c = v_g



% Gain and Detuning
B_env = 0 + 0*(1j);
B_gt = B_env.*ones(1, Ns);


% Grating
Kappa0 = zeros(1, Ns);
% Kappa_start = [1/8, 5/8];
% Kappa_end = [3/8, 7/8];

Kappa_start = [2/8];
Kappa_end = [6/8];
% k_B = 10000;
% x_B = k_B/(2*pi);
x_B = 40000;
k_B = x_B*2*pi;
d_er = 0.1;

for ii = (1:(size(Kappa_start)))
    if Kappa_end(ii) - Kappa_start(ii) > 0
        kp_s_ind = round(Kappa_start(ii)*Ns);
        kp_end_ind = round(Kappa_end(ii)*Ns);
        Kappa0(kp_s_ind:kp_end_ind) = kc*d_er/4;
    end
end


Kf_detune = exp(-(1j)*4*pi*(x_B - xc).*z_domain);
Kr_detune = exp((1j)*4*pi*(x_B - xc).*z_domain);


% Gain Dispersion
g_fwhm = 3.35*10^12;
Lgamma = 2*pi*g_fwhm;
Lw0 = 0;
Lgain = 0;
Cw0 = -Lgamma + 1j*Lw0;




% Gaussian params, defined in terms of ds
u_gaus = 250;
s_gaus = 50;
amp_gaus = 20;

% Input stream
input_length = 2*u_gaus + s_gaus*0 + 1;

amp_trig = 20;
freq = 100;
sample_rate = 1/dt;
% input_trig = f_trig(amp_trig, freq)


u_gaus_z = u_gaus * dz;
s_gaus_z = s_gaus * dz;

input_domain = linspace(0, input_length*dz, input_length);

E_in_static = zeros(1, input_length);
E_in_static(1, 1:input_length) = f_gaus(amp_gaus, u_gaus_z, s_gaus_z, input_domain);

E_in = E_in_static(1, 1:input_length);

mod_k = 50000;


E_in_mod = zeros(1, input_length);
E_in_mod(1, 1:input_length) = f_trig(amp_gaus, mod_k, input_domain);

% E_in = E_in.*real(E_in_mod(1, 1:input_length));



MAX_IN_VAL = max(E_in);


y_max = 1.2*max(E_in);

% Now we create the figure that will hold the graphs
fig_graph = figure("Name","Graphs");

% ===== (1) Top Left Plot - Holds the forward wave component =====
ax = subplot(3, 2, 1);
plt1 = plot(z_domain, real(Ef(1, :)), "-", z_domain, imag(Ef(1, :)), "--");
title("Forward Component");
xlabel("z (m)");
ylabel("E_f(z)");
xticks(dz*[200, 400, 600, 800, 1000]);
ylim([-y_max/2, y_max]);

% =====


% ===== (3) Middle Left Plot - Holds the reverse wave component =====
subplot(3, 2, 3);
plt3 = plot(z_domain, real(Er(1, :)), "-", z_domain, imag(Er(1, :)), "--");
title("Reverse Component");
xlabel("z (m)");
ylabel("E_r(z)");
xticks(dz*[200, 400, 600, 800, 1000]);
ylim([-y_max/2, y_max]);



% ===== (5) Bottom Left Plot - Wave Amplitudes at the boundaries =====
subplot(3, 2, 5);

left_input = zeros(1, Nt);
left_output = zeros(1, Nt);
right_input = zeros(1, Nt);
right_output = zeros(1, Nt);

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
fft_len = 2000;
% fft_len = double(Nt - 1);
freq_precision = sample_rate/fft_len;
total_time = fft_len/sample_rate;

% Frequencies represented in the fft
fft_domain = freq_precision*(-fft_len/2:(fft_len/2)-1);


% w_dom = 2.0*pi*(0:(t_sim/dt)-1)'/t_sim;
% kv = find(w_dom >= pi/dt);
% w_dom(kv) = w_dom(kv) - 2*pi/dt;

% Select how much of the fft to graph, basically selecting the "zoom",
% centered at the zero frequency
view_hw = 0.1;
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
fft_in_abs = abs(E_fft_in);

% These are what actually get plotted
fft_plot1 = fft_left_abs;
fft_plot2 = fft_right_abs;

fft_plot_in = fft_in_abs;

plt4 = plot( ...
    fft_domain(fft_lb:fft_ub), fft_plot1(fft_lb:fft_ub ), "-r", ...
    fft_domain(fft_lb:fft_ub), fft_plot2(fft_lb:fft_ub ), "-b", ...
    fft_domain(fft_lb:fft_ub), fft_plot_in(fft_lb:fft_ub), "--k");

title("Frequency Spectrum");
xlabel("Frequency (Hz)");
ylabel("E");

fft_min = min(min(abs(fft_plot1), abs(fft_plot2)));
fft_max = max(max(abs(fft_plot1), abs(fft_plot2)));

E_in_fft = fftshift(fft(E_in));
fft_min = min(abs(E_in_fft));
fft_max = max(abs(E_in_fft));

if fft_min == fft_max
    fft_min = fft_min - 1;
    fft_max = fft_max + 1;
end

% ylim([min(fft_in_abs), max(fft_in_abs)]);
ylim([fft_min, fft_max]);



% ===== (6) Bottom Right Plot - Phase Plot =====
subplot(3, 2, 6);

phase_plot_in = angle(E_fft_in);
phase_plot1 = angle(E_fft_left);
phase_plot2 = angle(E_fft_right);

angles1 = angle(E_fft_left);
angles2 = angle(E_fft_right);

total_phase2 = zeros(1, fft_len);
phase2_sum = 0;

angles_diff_1 = angles1(1, 2:fft_len) - angles1(1, 1:fft_len-1);
angles_diff_2 = angles2(1, 2:fft_len) - angles2(1, 1:fft_len-1);

% phase_plot1

% Select how much of the fft to graph, basically selecting the "zoom",
% centered at the zero frequency
view_phw = 0.4;
phase_lb = int32((0.5-view_phw)*fft_len);
phase_ub = int32((0.5+view_phw)*fft_len);

phase_plot_in_pos = zeros(1, fft_len);
phase_plot_in_pos(1, 1001:2000) = phase_plot_in(1, 1001:2000) + ...
flip(phase_plot_in(1, 1:1000));

phase_plot1_pos = zeros(1, fft_len);
phase_plot1_pos(1, 1001:2000) = phase_plot1(1, 1001:2000) + ...
flip(phase_plot1(1, 1:1000));

phase_plot2_pos = zeros(1, fft_len);
phase_plot2_pos(1, 1001:2000) = phase_plot2(1, 1001:2000) + ...
flip(phase_plot2(1, 1:1000));

phase1_diff = zeros(1, fft_len);
phase2_diff = zeros(1, fft_len);

plt6 = plot( ...
    fft_domain(phase_lb:phase_ub), phase_plot1(phase_lb:phase_ub), "-r", ...
    fft_domain(phase_lb:phase_ub), phase_plot2(phase_lb:phase_ub), "-b", ...
    fft_domain(phase_lb:phase_ub), phase_plot_in(phase_lb:phase_ub), "--k");


title("Phase Shift");
xlabel("Frequency (Hz)");
ylabel("Phase (E)");



Ef_test = zeros(1, Nz);
Er_test = zeros(1, Nz);


u_gaus_t = u_gaus * dt;
s_gaus_t = s_gaus * dt;

fft_len_z = 1000;
fft_precision_z = double(1/z_len);


num_chunks = size(Kappa_start, 1);

Ef_fft_chunks = zeros(num_chunks, fft_len_z);
Er_fft_chunks = zeros(num_chunks, fft_len_z);


% Frequencies represented in the spatial fft
fft_domain_z = fft_precision_z*(-fft_len_z/2:(fft_len_z/2)-1);

Ef_fft_z = zeros(1, int32(fft_len_z));
Er_fft_z = zeros(1, int32(fft_len_z));


% ========== SECOND GRAPH ==========
% ==================================

% Now we create the figure that will hold the graphs
fig2 = figure("Name","Graphs 2");
ax2 = subplot(2, 2, 1);

% ===== (7) Top Left Plot - Holds the forward wave component =====
plt7 = plot(z_domain, real(Ef(1, :)), "-", z_domain, imag(Ef(1, :)), "--");
title("Forward Component");
xlabel("z (m)");
ylabel("E_f(z)");
xticks(dz*[200, 400, 600, 800, 1000]);
ylim([-y_max, y_max]);



% ===== (9) Bottom Left Plot - Holds the reverse wave component =====
subplot(2, 2, 3);
plt9 = plot(z_domain, real(Er(1, :)), "-", z_domain, imag(Er(1, :)), "--");
% plt9 = plot(fft_domain_z, real(Ef_fft_z2(1, :)), "-", fft_domain_z, imag(Ef_fft_z2(1, :)), "--");
title("Reverse Component");
xlabel("z (m)");
ylabel("E_r(z)");
xticks(dz*[200, 400, 600, 800, 1000]);
ylim([-y_max, y_max]);


% ===== (8) Top Right Plot - Holds the reassembled forward wave component =====
subplot(2, 2, 2);
plt8 = plot(z_domain, real(Ef_test(1, :)), "-", z_domain, imag(Ef_test(1, :)), "--");
title("New Forward Component");
xlabel("z (m)");
ylabel("E_f(z)");
xticks(dz*[200, 400, 600, 800, 1000]);
ylim([-y_max, y_max]);



% ===== (10) Bottom Right Plot - Holds the reassembled reverse wave component =====
subplot(2, 2, 4);
plt10 = plot(z_domain, real(Er_test(1, :)), "-", z_domain, imag(Er_test(1, :)), "--");
% plt10 = plot(fft_domain_z, real(Ef_fft_z(1, :)), "-", fft_domain_z, imag(Ef_fft_z(1, :)), "--");
title("New Reverse Component");
xlabel("z (m)");
ylabel("E_r(z)");
xticks(dz*[200, 400, 600, 800, 1000]);
ylim([-y_max, y_max]);




pause_freq = 40;
break_freq = 10;
check_limits_freq = 5;

update_cmp_graph_freq = 3;


% Main simulation
for iter = 1:Nt
    
    % ===== (1) Step forward field ===== %

    % Save the last part of the segment for the reverse field
    last_fwd_seg = Ef(1, Ns);

    % Iterate the forward field component
    Ef(1, 2:Ns) = Ef(1, 1:Ns-1);

    % Grating
    Ef(1, 1:Ns-1) = Ef(1, 1:Ns-1) + ...
    Er(1, 2:Ns).*(1j).*Kappa0(2:Ns).*dz;


    % Ef_test(1, 1:Ns-1) = Ef_test(1, 1:Ns-1) + ...
    % Er_test(1, 2:Ns).*(1j).*Kappa0(2:Ns).*dz;

    % Gain and detuning
    Ef(1, 2:Ns) = Ef(1, 2:Ns) - ...
    Ef(1, 2:Ns).*(1j).*B_gt(1, 2:Ns).*dz;

    Pf(1, 2:Ns) = ( Pfp(1, 2:Ns).*(1+0.5*dt*Cw0) + ...
0.5*dt*Lgamma.*(Ef(1, 2:Ns)+Efp(1, 2:Ns)) )./(1-0.5*dt*Cw0);

    Ef(1, 2:Ns) = Ef(1, 2:Ns) - ...
    Lgain*(Ef(1, 2:Ns) - Pf(1, 2:Ns));


    % Update the left side of the forward travelling wave
    Ef(1, 1) = E_in(1, input_length);


    set(plt1(1), 'XData', z_domain, 'YData', real(Ef));
    set(plt1(2), 'XData', z_domain, 'YData', imag(Ef));




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

    fft_in_abs = abs(E_fft_in);


    fft_plot_in = E_fft_in;
    fft_plot1 = E_fft_left;
    fft_plot2 = E_fft_right;

    set(plt4(1), 'XData', fft_domain(fft_lb:fft_ub), 'YData', abs(fft_plot1(fft_lb:fft_ub)));
    set(plt4(2), 'XData', fft_domain(fft_lb:fft_ub), 'YData', abs(fft_plot2(fft_lb:fft_ub)));
    set(plt4(3), 'XData', fft_domain(fft_lb:fft_ub), 'YData', abs(fft_plot_in(fft_lb:fft_ub)));


    % Update phase
    phase_plot1 = angle(E_fft_left);
    phase_plot2 = angle(E_fft_right);

    angles1 = angle(E_fft_left);
    angles2 = angle(E_fft_right);

    angles2_negs = find(angles2 < 0);
    angles2(angles2_negs) = angles2(angles2_negs) + 2*pi;

    % phase2_sum = phase2_sum + angles2()
    % total_phase2(1, iter) = 
    
    angles_diff_1 = angles1(1, 2:fft_len) - angles1(1, 1:fft_len-1);
    angles_diff_2 = angles2(1, 2:fft_len) - angles2(1, 1:fft_len-1);

    phase_plot_in = angle(E_fft_in);

    phase1_diff = phase_plot1 - phase_plot_in;
    phase2_diff = phase_plot2 - phase_plot_in;

    % set(plt6(1), 'XData', fft_domain(phase_lb:phase_ub), 'YData', phase1_diff(phase_lb:phase_ub));
    % set(plt6(2), 'XData', fft_domain(phase_lb:phase_ub), 'YData', phase2_diff(phase_lb:phase_ub));
    set(plt6(3), 'XData', fft_domain(phase_lb:phase_ub), 'YData', phase_plot_in(phase_lb:phase_ub));

    set(plt6(1), 'XData', fft_domain(phase_lb:phase_ub), 'YData', phase_plot1_pos(phase_lb:phase_ub));
    set(plt6(2), 'XData', fft_domain(phase_lb:phase_ub), 'YData', phase_plot2_pos(phase_lb:phase_ub));


    % ===== Update reconstructed/actual comparison graphs =====
    if mod(iter, update_cmp_graph_freq) == 0
        % set(plt7, 'XData', fft_domain(phase_lb:phase_ub), 'YData', phase_plot1(phase_lb:phase_ub)*0);

    end


    % ===== (3) Step reverse field =====
    Er(1, 1:Ns-1) = Er(1, 2:Ns);
    
    % Gain and detuning
    Er(1, 1:Ns-1) = Er(1, 1:Ns-1) - ...
    Er(1, 1:Ns-1).*(1j).*B_gt(1, 1:Ns-1).*dz;

    % Grating
    Er(1, 2:Ns) = Er(1, 2:Ns) + ...
    Ef(1, 1:Ns-1).*(1j).*Kappa0(1:Ns-1).*dz;


    % Dispersion
    Pr(1, 2:Ns) = ( Prp(1, 2:Ns)*(1+0.5*dt*Cw0) + ...
0.5*dt*Lgamma*(Er(1, 2:Ns)+Erp(1, 2:Ns)) )./(1-0.5*dt*Cw0);

    Er(1, 2:Ns) = Er(1, 2:Ns) - ...
    Lgain*(Er(1, 2:Ns) - Pr(1, 2:Ns));


    set(plt3(1), 'XData', z_domain, 'YData', real(Er(1, :)));
    set(plt3(2), 'XData', z_domain, 'YData', imag(Er(1, :)));


    % Step new forward field
    gr = Kappa0(1, :) .* Kr_detune(1, :) .* Ef_test(1, :);
    % Ef_test(1, 2:Nz) = Ef_test(1, 1:(Nz-1));


    end_index = Nz-1;
    Ef_fft_z = fftshift(fft(Ef_test(1, 1:end_index)));
    Er_fft_z = fftshift(fft(Er_test(1, 2:end_index+1)));

    Ef_curr_chunk = zeros(1, Nz);
    Er_curr_chunk = zeros(1, Nz);

    gf = Kappa0 .* Kf_detune .* Er_test;

    Gf_fft_z = fftshift(fft(gf(1, 1:end_index)));
    Gr_fft_z = fftshift(fft(gr(1, 1:end_index)));

    Ef_fft_z = Ef_fft_z + (1j)*dz*Gf_fft_z;
    Ef_test(1, 2:Nz) = ifft(fftshift(Ef_fft_z));

    Er_fft_z = Er_fft_z + (1j)*dz*Gr_fft_z;
    Er_test(1, 1:(Nz-1)) = ifft(fftshift(Er_fft_z));


%{
    Energyf = sum(abs(Ef_curr_chunk).^2)/fft_len_z;
    Energyr = sum(abs(Er_curr_chunk).^2)/fft_len_z;

    ESDf = (abs(Ef_curr_chunk).^2)/(Energyf*Nz);
    ESDr = (abs(Er_curr_chunk).^2)/(Energyr*Nz);

    esdf_hbw = 0.05;
    ESDf_lb = int32(fft_len_z*(0.5-esdf_hbw));
    ESDf_ub = int32(fft_len_z*(0.5+esdf_hbw));
    ESDf_inner = ESDf(ESDf_lb:ESDf_ub);

    inner_energy_frac = sum(ESDf_inner);
    inner_energy_threshold = 0.97;

    while (inner_energy_frac < inner_energy_threshold)
        esdf_hbw = esdf_hbw + 0.05;

        ESDf_lb = int32(fft_len_z*(0.5-esdf_hbw));
        ESDf_ub = int32(fft_len_z*(0.5+esdf_hbw));
        ESDf_inner = ESDf(ESDf_lb:ESDf_ub);

        inner_energy_frac = sum(ESDf_inner);
    end
    %}



    Ef_test(1, 1) = E_in(1, input_length);

    % (4) Finally, step the input function
    E_in(1, 2:input_length) = E_in(1, 1:(input_length-1));
    E_in(1, 1) = 0;


    % (5) Update the previous fields
    Efp = Ef;
    Erp = Er;
    Pfp = Pf;
    Prp = Pr;



    % Update second figure graphs
    set(plt7(1), 'XData', z_domain, 'YData', real(Ef));
    set(plt7(2), 'XData', z_domain, 'YData', imag(Ef));


    set(plt9(1), 'XData', z_domain, 'YData', real(Er));
    set(plt9(2), 'XData', z_domain, 'YData', imag(Er));


    set(plt8(1), 'XData', z_domain, 'YData', real(Ef_test));
    set(plt8(2), 'XData', z_domain, 'YData', imag(Ef_test));

    set(plt10(1), 'XData', z_domain, 'YData', real(Er_test));
    set(plt10(2), 'XData', z_domain, 'YData', imag(Er_test));


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
    trig_func = amp.*exp(-(1j)*2*pi*freq*t);
end

function approx_g = f_gaus2(ug, sg, f)
    % p1 = -2*pi*(vp0^2).*f;
    % p2 = pi*(sg^2).*f - (1j)*vp0.*(ug.*f + z);
    % temp1 = (ug.*f + z);
    % test44 = ((1j).*vp0).*(temp1);
    % approx_g = (1./vp0)*exp(-(2*pi*(vp0^2)).*f .* ((pi*(sg^2).*f - ((1j)*vp0).*(ug.*f + z) )) );
    gaus_mag = 1/(sqrt(2*pi)*sg);
    approx_g = gaus_mag*exp(-(2*(pi)*f).*(pi*(sg^2)*f + ((1j)*ug) ) );
end

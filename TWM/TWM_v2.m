
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

vg = vp0;

% Carrier Characteristics
lambda_c = 25*dz;
xc = 1/lambda_c;
kc = 2*pi*xc;
fc = vp0/lambda_c;
wc = 2*pi*fc;


% Gain and Detuning
B_env = 0 + 0*(1j);
B_gt = B_env.*ones(1, Ns);


% Grating
Kappa0 = zeros(1, Ns);

Kappa_start = [2/16, 8/16];
Kappa_end = [4/16, 12,16];

% Bragg Wavelength
lambda_b = 22*dz;
xb = 1/lambda_b;
kb = xb*2*pi;
fb = vp0/lambda_b;
wb = 2*pi*fb;

d_er = 0.4;

Kappa_mag = (kc*d_er)/4;

for ii = (1:(size(Kappa_start)))
    if Kappa_end(ii) - Kappa_start(ii) > 0
        kp_s_ind = round(Kappa_start(ii)*Ns);
        kp_end_ind = round(Kappa_end(ii)*Ns);
        Kappa0(kp_s_ind:kp_end_ind) = Kappa_mag;
    end
end


Kf_detune = exp(-(1j)*4*pi*(xb - xc).*z_domain);
Kr_detune = exp((1j)*4*pi*(xb - xc).*z_domain);


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
freq_precision = sample_rate/fft_len;
total_time = fft_len/sample_rate;

% Frequencies represented in the fft
fft_domain = freq_precision*(-fft_len/2:(fft_len/2)-1);


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


% Select how much of the fft to graph, basically selecting the "zoom",
% centered at the zero frequency
view_phw = 0.4;
phase_lb = int32((0.5-view_phw)*fft_len);
phase_ub = int32((0.5+view_phw)*fft_len);


plt6 = plot( ...
    fft_domain(phase_lb:phase_ub), phase_plot1(phase_lb:phase_ub), "-r", ...
    fft_domain(phase_lb:phase_ub), phase_plot2(phase_lb:phase_ub), "-b", ...
    fft_domain(phase_lb:phase_ub), phase_plot_in(phase_lb:phase_ub), "--k");


title("Phase Shift");
xlabel("Frequency (Hz)");
ylabel("Phase (E)");


u_gaus_t = u_gaus * dt;
s_gaus_t = s_gaus * dt;

fft_len_z = 1000;
fft_precision_z = double(1/z_len);


% Frequencies represented in the spatial fft
fft_domain_z = fft_precision_z*(-fft_len_z/2:(fft_len_z/2)-1);

Ef_fft_z = zeros(1, int32(fft_len_z));
Er_fft_z = zeros(1, int32(fft_len_z));


pause_freq = 40;
break_freq = 10;
check_limits_freq = 5;

update_cmp_graph_freq = 3;


% Conservation of energy term
fE = sqrt(1-(Kappa_mag*dz)^2);

vg_cf = ones(1, Ns);

cf_iters = 3;



% Main simulation
for iter = 1:Nt
    
    for kk = 1:cf_iters
    fc = vg/lambda_c;
    wc = 2*pi*fc;
    
    fb = vg/lambda_b;
    wb = 2*pi*fb;
    
    w_diff = wc - wb;
    if wc == wb
        vg_cf = ones(1, Ns);
    else
        vg_cf = sqrt((w_diff^2 - Kappa0.^2)/(w_diff^2));
    end

    vg = vp0.*vg_cf;
    end


    % Update vg


    % ===== (1) Step forward field ===== %

    % Save the last part of the segment for the reverse field
    last_fwd_seg = Ef(1, Ns);

    %{
    % Iterate the forward field component
    Ef(1, 2:Ns) = Ef(1, 1:Ns-1);

    % Grating
    Ef(1, 1:Ns-1) = Ef(1, 1:Ns-1) + ...
    Er(1, 2:Ns).*(1j).*Kappa0(2:Ns).*dz;

    %}

    end_index = Nz-1;

    % Update forward field
    Ef = fE*Ef;


    % Gain and detuning
    Ef(1, 2:Ns) = Ef(1, 2:Ns) - ...
    Ef(1, 2:Ns).*(1j).*B_gt(1, 2:Ns).*dz;
    
    % Grating
    Ef_fft_z = fftshift(fft(Ef(1, 1:end_index)));
    gf = Kappa0 .* Kf_detune .* Er;
    Gf_fft_z = fftshift(fft(gf(1, 1:end_index)));
    Ef_fft_z = Ef_fft_z + (1j)*dz*Gf_fft_z;
    Ef(1, 2:Nz) = ifft(fftshift(Ef_fft_z));

    % Dispersion

    % Update the Pf field using the updated Ef field
    Pf(1, 2:Ns) = ( Pfp(1, 2:Ns).*(1+0.5*dt*Cw0) + ...
0.5*dt*Lgamma.*(Ef(1, 2:Ns)+Efp(1, 2:Ns)) )./(1-0.5*dt*Cw0);

    % 
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

    phase_plot_in = angle(E_fft_in);

    set(plt6(3), 'XData', fft_domain(phase_lb:phase_ub), 'YData', phase_plot_in(phase_lb:phase_ub));

    set(plt6(1), 'XData', fft_domain(phase_lb:phase_ub), 'YData', phase_plot1(phase_lb:phase_ub));
    set(plt6(2), 'XData', fft_domain(phase_lb:phase_ub), 'YData', phase_plot2(phase_lb:phase_ub));


    % ===== Update reconstructed/actual comparison graphs =====
    if mod(iter, update_cmp_graph_freq) == 0
        % set(plt7, 'XData', fft_domain(phase_lb:phase_ub), 'YData', phase_plot1(phase_lb:phase_ub)*0);

    end


    % ===== (3) Step reverse field =====


    % Gain and detuning
    Er(1, 1:Ns-1) = Er(1, 1:Ns-1) - ...
    Er(1, 1:Ns-1).*(1j).*B_gt(1, 1:Ns-1).*dz;

    
    % Bragg Grating
    Er = fE*Er;
    Er_fft_z = fftshift(fft(Er(1, 2:end_index+1)));
    gr = Kappa0(1, :) .* Kr_detune(1, :) .* Ef(1, :);
    Gr_fft_z = fftshift(fft(gr(1, 1:end_index)));
    Er_fft_z = Er_fft_z + (1j)*dz*Gr_fft_z;
    Er(1, 1:(Nz-1)) = ifft(fftshift(Er_fft_z));


    % Dispersion

    % Update the Pr field using the updated Er field
    Pr(1, 2:Ns) = ( Prp(1, 2:Ns)*(1+0.5*dt*Cw0) + ...
0.5*dt*Lgamma*(Er(1, 2:Ns)+Erp(1, 2:Ns)) )./(1-0.5*dt*Cw0);

    Er(1, 2:Ns) = Er(1, 2:Ns) - ...
    Lgain*(Er(1, 2:Ns) - Pr(1, 2:Ns));


    set(plt3(1), 'XData', z_domain, 'YData', real(Er(1, :)));
    set(plt3(2), 'XData', z_domain, 'YData', imag(Er(1, :)));



    % (4) Finally, step the input function
    E_in(1, 2:input_length) = E_in(1, 1:(input_length-1));
    E_in(1, 1) = 0;


    % (5) Update the previous fields
    Efp = Ef;
    Erp = Er;
    Pfp = Pf;
    Prp = Pr;


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

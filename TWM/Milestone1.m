% Define start and end points (Say these are in meters)
X1 = 0;
X2 = 2;
res = 0.025;
set_length = ((X2-X1)/res) + 1;

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


% Input stream
input_length = 10;
E_in = zeros(1, input_length);
E_in(1, 1:input_length) = f_gaus(0:input_length-1);



% Now we create the figure that will hold the graphs
fig_graph = figure("Name","Graphs");

% Top Plot - Holds the forward wave component
ax = subplot(3, 1, 1);
plt1 = plot(x, real(Ef(1, :)), "-", x, imag(Ef(1, :)), "--");
title("Forward Component");
xlabel("z (m)");
ylabel("E_f(z)");


% Middle Plot - Holds the forward wave component
subplot(3, 1, 2);
plt2 = plot(x, real(Er(1, :)), "-", x, imag(Er(1, :)), "--");
title("Reverse Component");
xlabel("z (m)");
ylabel("E_r(z)");



% Bottom Plot - Wave Input/Output Amplitudes at the boundaries
subplot(3, 1, 3);

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

legend("Left Input", "Left Output", "Right Input", "Right Output");



% Main simulation
for iter = 1:sim_ticks
    
    % ===== (1) Step forward field ===== %

    % Save the last part of the segment for the reverse field
    last_fwd_seg = Ef(1, set_length);

    % Iterate the forward field component
    Ef(1, 2:set_length) = Ef(1, 1:set_length-1);

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



    % ===== (3) Step reverse field =====
    Er(1, 1:set_length-1) = Er(1, 2:set_length);
    Er(1, set_length) = last_fwd_seg * RR;

    set(plt2(1), 'XData', x, 'YData', real(Er(1, :)));
    set(plt2(2), 'XData', x, 'YData', imag(Er(1, :)));



    % (4) Finally, step the input function
    E_in(1, 2:input_length) = E_in(1, 1:(input_length-1));
    E_in(1, 1) = 0;
    

    pause(0.001);

end



% Generate gaussian pdf
function func_gaussian = f_gaus(t)
    u_gaus = 5;
    s_gaus = 2;
    func_gaussian = (1./(s_gaus.*sqrt(2.*pi))).*exp( -0.5.*((t - u_gaus)./s_gaus).^2 );
end


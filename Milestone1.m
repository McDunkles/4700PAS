% Define start and end points (Say these are in meters)
X1 = 0;
X2 = 10;
res = 0.1;
set_length = ((X2-X1)/res) + 1;


% Wave propagation speed
c = 3 * 10^8;
n = 3;
vg = c/n;


% Simulation step time, in seconds
dt = 0.5 * (10^-9);

% dz-dt Synchronization, dz given in meters
dz = vg*dt;
frac_z = 0;


% Reflection Coefficients
RL = 0.9*(1j);
RR = 0.9*(1j);


% 1D spatial domain
x = linspace(X1,X2,set_length);

% Starting time is 0
t = 0;

% Forward and reverse wave components
Ef = zeros(1, set_length);
Er = zeros(1, set_length);


% Input stream
input_length = 10;
E_in = zeros(1, input_length);


% Populate the values to feed in as input
for iter = 1:input_length
    E_in(1, iter) = iter*(10-iter);
end


% Create the control panel
fig_width = 560;
fig_height = 420;

global fig
fig = uifigure("Name","Control Panel");
fig.Position = [fig.Position(1:2), fig_width, fig_height];
fig.Visible = "off";

% Toggle the control panel's visibility by pressing 'v' while the
% simulation is running
% (toggle_vis() is defined at the bottom of this script)
set(fig, "KeyPressFcn", @toggle_vis);

% Root layout, which will hold the 3 sublayouts (header, body, and footer)
gl_root = uigridlayout(fig);
gl_root.RowHeight = {"1x", "3x", "1x"};
gl_root.ColumnWidth = {"1x"};
gl_root.ColumnSpacing = 0;
gl_root.RowSpacing = 0;
gl_root.Padding = 0; % [0, 0, 0, 0]

% Header layout
gl_header = uigridlayout(gl_root);
gl_header.RowHeight = {"1x"};
gl_header.ColumnWidth = {"1x"};
gl_header.Layout.Row = 1;
gl_header.BackgroundColor = "red";

% Body layout
gl_body = uigridlayout(gl_root);
gl_body.RowHeight = {"1x", "2x"};
gl_body.ColumnWidth = {"1x", "1x", "1x"};
gl_body.Layout.Row = 2;
gl_body.BackgroundColor = "blue";

% Footer layout
gl_footer = uigridlayout(gl_root);
gl_footer.RowHeight = {"1x"};
gl_footer.ColumnWidth = {"1x"};
gl_footer.Layout.Row = 3;
gl_footer.BackgroundColor = "green";

% Title label
lb_title = uilabel(gl_header);
lb_title.Text = "Control Panel";
lb_title.HorizontalAlignment = "center";
lb_title.FontSize = 30;

% Sadly this does not seem to work, I may look into this, but
% it's not a priority
lb_title.FontName = "Comic Sans";

% Slider to control the index of refraction (n) of the medium,
% pretty much just a speed modifier right now
sld = uislider(gl_body);
sld.Layout.Row = 1;
sld.Layout.Column = 2;

% Initialize the slider value and other properties of the slider
sld.Value = n;
n_bounds = [1, 6];
set(sld, "Limits", n_bounds);
sld.MajorTicks = 1:6;
sld.MinorTicks = (1:12).*0.5;

% Label for the slider
lb = uilabel(gl_body);
lb.Text = "Index of Refraction (n)";
lb.Layout.Row = 1;
lb.Layout.Column = 1;
lb.FontSize = 16;



% Now we create the figure that will hold the graphs
fig_graph = figure("Name","Graphs");

% Bind the 'toggle_vis' function to this figure too so we can
% make the control panel visible again
set(fig_graph, "KeyPressFcn", @toggle_vis);












ax = subplot(3, 1, 1);
% properties(ax)


% hold on
size_ef_real = size(real(Ef(1, :)))
size_ef_imag = size(imag(Ef(1, :)))
size_er_real = size(real(Er(1, :)))
size_er_imag = size(imag(Er(1, :)))
plt1 = plot(x, real(Ef(1, :)), "-", x, imag(Ef(1, :)), "--");

subplot(3, 1, 2);
plt2 = plot(x, real(Er(1, :)), "-", x, imag(Er(1, :)), "--");



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
for iter = 1:10*num_steps

    n = sld.Value;
    vg = c/n;
    dz = vg*dt;

    % Compute shift and fraction values
    shift_fl = dz/res;
    shift = cast(shift_fl, "int32");
    frac_z = frac_z + (shift_fl - double(shift));
    
    if (frac_z >= 0.5)
        shift = shift + 1;
        frac_z = frac_z - 1;

    elseif (frac_z <= -0.5)
            shift = shift - 1;
            frac_z = frac_z + 1;
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
    if (iter < 20)

        msg_format = "[%d] Shift = %d, frac_z = %f, total_shift = %d, prev_total_shift = %d\n";
        message = sprintf(msg_format, iter, shift, frac_z, total_shift, prev_total_shift);

        fprintf(message);
    end

    
    % ===== Step forward field ===== %
    last_fwd_seg = Ef(1, shift_comp:set_length);

    % Translate all elements over by 'shift' units
    Ef(1, 1+shift:set_length) = Ef(1, 1:set_length-shift);

    % Zero out the outlying first 'shift' rows, ready to be written with
    % new data

    Er_le = flip(Er(1, 1:shift), 2);

    Ef_real_le = real(Er_le);
    Ef_imag_le = imag(Er_le);
    Ef(1, 1:shift) = Ef_real_le(1, 1:shift) + (1j)*Ef_imag_le(1, 1:shift);

    Ef(1, 1:shift) = Ef(1, 1:shift) * RL;




    % If the input function is still going, add it to the beginning of Ef
    % Should be filling the first 'shift' elements
    if total_shift <= input_length
        Ef(1, 1:shift) = Ef(1, 1:shift) + flip(E_in(1, prev_total_shift+1:total_shift), 2);
    end


    Ef_real = real(Ef(1, :));
    Ef_imag = imag(Ef(1, :));

    

    set(plt1(1), 'XData', x, 'YData', Ef_real);
    set(plt1(2), 'XData', x, 'YData', Ef_imag);


    % Step reverse field
    Er(1, 1:set_length-shift) = Er(1, 1+shift:set_length);
    Er(1, shift_comp:set_length) = flip(last_fwd_seg(1, :), 2);

    Er(1, shift_comp:set_length) = Er(1, shift_comp:set_length) * RR;


    set(plt2(1), 'XData', x, 'YData', real(Er(1, :)));
    set(plt2(2), 'XData', x, 'YData', imag(Er(1, :)));
    

    pause(0.02);

prev_total_shift = total_shift;
end

global vis
vis = 1;

% Definition of the function for toggling the control panel visibility
% Very crudely written right now, will refactor later
function toggle_vis(src, event)
    global vis fig

    % See 'Constants.m' for the definition of the Constants class
    C = Constants;

    key = uint8(event.Key);
    if (key == C.V_ASCII)

        vis_sig = "";

        if vis == 1
            vis = 0;
            vis_sig = "off";
        else
            vis = 1;
            vis_sig = "on";
        end

        set(fig, "Visible", vis_sig);

    end


end

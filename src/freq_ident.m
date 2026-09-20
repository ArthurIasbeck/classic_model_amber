clc; clear; close all; 

data = load("../data/open_loop_frequency_responses.txt");

w = data(:, 1);
Gjw = data(:, 2) + 1j * data(:, 3);
data = idfrd(Gjw, w, 0);

sys = tfest(data, 5);

opts = bodeoptions;
opts.PhaseWrapping = 'on';
opts.PhaseWrappingBranch = 0;
bode(data, sys, w, opts)
[mag_data, phase_data, w_data] = bode(data, w);
[mag_sys, phase_sys, ~] = bode(sys, w);

mag_data = squeeze(mag_data);
phase_data = squeeze(phase_data);
mag_sys = squeeze(mag_sys);
phase_sys = squeeze(phase_sys);
w_data = squeeze(w_data);
writematrix([w_data(:), mag_data(:), phase_data(:), mag_sys(:), phase_sys(:)], ...
    "../data/ident_freq_v13_freq.txt");

csiN1_W13_DE = 0.08;
fN1_W13_DE = 2.324778563656447e3;
csiD1_W13_DE = 0.37;
fD1_W13_DE = 1.759291886010284e3;
csiN2_W13_DE = 0.03;
fN2_W13_DE = 4.712388980384690e3;
csiD2_W13_DE = 0.16;
fD2_W13_DE = 4.178318229274425e3;

s = tf('s');

C = 0.0062*(1 + 46/s)*(400/77)*((s + 77*2*pi)/(s + 400*2*pi))*(409/124)*((s + 124*2*pi)/(s + 409*2*pi))*tf([1 2*csiN1_W13_DE*fN1_W13_DE fN1_W13_DE^2],[1 2*csiD1_W13_DE*fD1_W13_DE fD1_W13_DE^2])*tf([1 2*csiN2_W13_DE*fN2_W13_DE fN2_W13_DE^2],[1 2*csiD2_W13_DE*fD2_W13_DE fD2_W13_DE^2])*0.4;

T = feedback(C * sys, 1); 
U = feedback(C, sys);

load("../data/prbs_v13.mat")
t = prbs_v13{:, 1};
i_real = prbs_v13{:, 2};
d = prbs_v13{:, 6};
y_real = prbs_v13{:, 10};

t = t(~isnan(d));
y_real = y_real(~isnan(d));
i_real = i_real(~isnan(d));
d = d(~isnan(d));

t = linspace(t(1), t(end), size(t, 1));
y = lsim(T, -d, t);
i = lsim(U, -d, t);
i_real_plot = i_real - mean(i_real);
writematrix([t(:), y(:), y_real(:), i(:), i_real_plot(:)], ...
    "../data/ident_freq_v13_time.txt");

figure;
plot(t, y); hold on;
plot(t, y_real);
grid on;
legend("Model", "Experiment")

figure;
plot(t, i); hold on;
plot(t, i_real_plot);
grid on;
legend("Model", "Experiment")

clc; clear; close all; 

data = load("../data/open_loop_frequency_responses.txt");

w = data(:, 1);
Gjw = data(:, 2) + 1j * data(:, 3);
data = idfrd(Gjw, w, 0);

sys = tfest(data, 8, 5);
display(sys)

sys_tf = tf(sys);
[numerator, denominator] = tfdata(sys_tf, "v");

model.numerator = numerator;
model.denominator = denominator;
model.sample_time = 0;
model.time_unit = "seconds";

model_file = "../data/freq_ident_model.json";
file_id = fopen(model_file, "w");
if file_id == -1
    error("Could not open model output file: %s", model_file);
end
fprintf(file_id, "%s\n", jsonencode(model, "PrettyPrint", true));
fclose(file_id);

opts = bodeoptions;
opts.PhaseWrapping = 'on';
opts.PhaseWrappingBranch = 0;
bode(data, sys, w, opts)

G_zpk = zpk(sys);
display(G_zpk)

figure;
pzmap(G_zpk)

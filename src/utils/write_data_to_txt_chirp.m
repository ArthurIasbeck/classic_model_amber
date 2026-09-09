clc; clear; close all;

load(strcat('../../data/', 'chirp_v13'));
load(strcat('../../data/', 'chirp_w13'));
load(strcat('../../data/', 'chirp_v24'));
load(strcat('../../data/', 'chirp_w24'));

i_v13 = [];
i_w13 = [];
i_v24 = [];
i_w24 = [];
d_v13 = [];
d_v24 = [];
d_w13 = [];
d_w24 = [];
y_v13 = [];
y_v24 = [];
y_w13 = [];
y_w24 = [];
t_end = [];

dataset_names = {'chirp_v13', 'chirp_w13', 'chirp_v24', 'chirp_w24'};

for i = 1:4
    data = eval(dataset_names{i});

    t_end(end + 1) = data.X.Data(end);

    i_v13 = [i_v13, data.Y(1).Data];
    i_w13 = [i_w13, data.Y(2).Data];
    i_v24 = [i_v24, data.Y(3).Data];
    i_w24 = [i_w24, data.Y(4).Data];

    d_v13 = [d_v13, data.Y(5).Data];
    d_v24 = [d_v24, data.Y(6).Data];
    d_w13 = [d_w13, data.Y(7).Data];
    d_w24 = [d_w24, data.Y(8).Data];

    y_v13 = [y_v13, data.Y(9).Data];
    y_v24 = [y_v24, data.Y(10).Data];
    y_w13 = [y_w13, data.Y(11).Data];
    y_w24 = [y_w24, data.Y(12).Data];
end

N = size(y_v13, 2);
t_total = sum(t_end);
t = linspace(0, t_total, N);

array_data = [t', i_v13', i_w13', i_v24', i_w24', d_v13', d_w13', ...
    d_v24', d_w24', y_v13', y_w13', y_v24', y_w24'];

writematrix(array_data, 'txt_data/chirp.txt')

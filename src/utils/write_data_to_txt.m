clc; clear; close all;

file_name = 'chirp_v13_4';
load(strcat('../../data/', file_name));
data = chirp_v13_4;

t = data.X.Data;

i_v13 = data.Y(1).Data;
i_w13 = data.Y(2).Data;
i_v24 = data.Y(3).Data;
i_w24 = data.Y(4).Data;

plot_data(t, i_v13, i_w13, i_v24, i_w24, 'i [A]')

d_v13 = data.Y(5).Data;
d_v24 = data.Y(6).Data;
d_w13 = data.Y(7).Data;
d_w24 = data.Y(8).Data;

plot_data(t, d_v13, d_w13, d_v24, d_w24, 'd [um]')

y_v13 = data.Y(9).Data;
y_v24 = data.Y(10).Data;
y_w13 = data.Y(11).Data;
y_w24 = data.Y(12).Data;

plot_data(t, y_v13, y_w13, y_v24, y_w24, 'y [um]')

array_data = [t', i_v13', i_w13', i_v24', i_w24', d_v13', d_w13', ...
    d_v24', d_w24', y_v13', y_w13', y_v24', y_w24'];

writematrix(array_data, strcat('txt_data/', file_name))

function plot_data(t, v13, w13, v24, w24, label)
    fig = figure;
    fig.Position = [100 100 1500 800];
    subplot(2, 2, 1)
    plot(t, v13)
    xlabel('t [s]')
    ylabel(label)
    title('v_{13}')

    subplot(2, 2, 2)
    plot(t, w13)
    xlabel('t [s]')
    ylabel(label)
    title('w_{13}')

    subplot(2, 2, 3)
    plot(t, v24)
    xlabel('t [s]')
    ylabel(label)
    title('v_{24}')

    subplot(2, 2, 4)
    plot(t, w24)
    xlabel('t [s]')
    ylabel(label)
    title('w_{24}')
end

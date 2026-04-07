%% plot_polefigures.m
% Построение полюсных фигур {111}, {110}, {100}
clear; clc; close all;

addpath('D:\MATLAB\MTEX\mtex-develop');
startup_mtex;

cs = crystalSymmetry('m-3m');
h = [Miller(1,1,1,cs), Miller(1,1,0,cs), Miller(1,0,0,cs)];

loadings = {'simple_shear', 'uniaxial_tension', 'uniaxial_compression', 'plane_strain'};
loading_labels = {'Simple Shear', 'Uniaxial Tension', 'Uniaxial Compression', 'Plane Strain'};
strains = [30, 50, 100];

% ======================== ОТДЕЛЬНЫЕ ФИГУРЫ ========================

for li = 1:length(loadings)
    loading = loadings{li};
    label   = loading_labels{li};

    for si = 1:length(strains)
        strain = strains(si);
        filename = sprintf('orientations_%s_%dpct.txt', loading, strain);

        if ~exist(filename, 'file')
            fprintf('SKIP: %s\n', filename);
            continue;
        end

        data = readmatrix(filename);
        qw = data(:, 1); qx = data(:, 2); qy = data(:, 3); qz = data(:, 4);

        ori = orientation(quaternion(qw, qx, qy, qz), cs);
        odf = calcDensity(ori, 'halfwidth', 5*degree);

        fprintf('%s %d%%: TI = %.2f\n', label, strain, textureindex(odf));

        fig = figure;
        plotPDF(odf, h, 'contourf', 'antipodal');
        mtexColorbar;
        title(sprintf('%s — %d%% strain', label, strain), 'FontSize', 16);

        outname = sprintf('PF_%s_%dpct.png', loading, strain);
        saveas(fig, outname);
        close(fig);
        fprintf('  Saved: %s\n', outname);
    end
end

% ======================== СРАВНИТЕЛЬНАЯ (100%) ========================

fprintf('\n=== Comparison 100%% ===\n');

fig = figure('Position', [100, 100, 1400, 1200]);

for li = 1:length(loadings)
    loading = loadings{li};
    label   = loading_labels{li};
    strain = 100;
    filename = sprintf('orientations_%s_%dpct.txt', loading, strain);

    data = readmatrix(filename);
    qw = data(:, 1); qx = data(:, 2); qy = data(:, 3); qz = data(:, 4);

    ori = orientation(quaternion(qw, qx, qy, qz), cs);
    odf = calcDensity(ori, 'halfwidth', 5*degree);

    subplot(4, 3, (li-1)*3 + 1);
    plotPDF(odf, h, 'contourf', 'antipodal');
    mtexColorbar;

    % Подпись строки
    subplot(4, 3, (li-1)*3 + 1);
    ylabel(label, 'FontSize', 12, 'Rotation', 0);
end

sgtitle(sprintf('Strain: %d%%', strain), 'FontSize', 16);
saveas(fig, sprintf('PF_compare_%dpct.png', strain));
close(fig);

fprintf('Done!\n');

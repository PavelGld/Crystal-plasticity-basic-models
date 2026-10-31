%% plot_texture.m -- pole figures for the UMAT implicit texture validation
clear; clc; close all;

addpath('D:\MATLAB\MTEX\mtex-develop');
startup_mtex;

cs = crystalSymmetry('m-3m');
h = [Miller(1,1,1,cs), Miller(1,1,0,cs), Miller(1,0,0,cs)];
root = 'D:/Pet_Projects/EigenProject/ver_PPP';

modes  = {'simple_shear', 'uniaxial_tension', 'uniaxial_compression', 'plane_strain'};
labels = {'Simple Shear', 'Uniaxial Tension', 'Uniaxial Compression', 'Plane Strain'};
pcts   = [30, 50, 100];

for i = 1:length(modes)
    for p = pcts
        f = fullfile(root, sprintf('texture_umat_%s_%dpct.txt', modes{i}, p));
        data = readmatrix(f, 'NumHeaderLines', 1);
        qw = data(:,2); qx = data(:,3); qy = data(:,4); qz = data(:,5);

        ori = orientation(quaternion(qw, qx, qy, qz), cs);
        odf = calcDensity(ori, 'halfwidth', 5*degree);

        fprintf('%s %d%%: TI = %.3f\n', labels{i}, p, textureindex(odf));

        fig = figure;
        plotPDF(odf, h, 'contourf', 'antipodal');
        mtexColorbar;
        title(sprintf('%s - %d%% (UMAT implicit)', labels{i}, p), 'FontSize', 14, 'Interpreter', 'none');
        outname = fullfile(root, sprintf('PF_umat_%s_%dpct.png', modes{i}, p));
        saveas(fig, outname);
        close(fig);
        fprintf('  saved %s\n', outname);
    end
end
fprintf('DONE\n');

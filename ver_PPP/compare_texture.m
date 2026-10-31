%% compare_texture.m -- explicit (existing validation) vs implicit (UMAT) textures
clear; clc; close all;
addpath('D:\MATLAB\MTEX\mtex-develop');
startup_mtex;

cs = crystalSymmetry('m-3m');
root = 'D:/Pet_Projects/EigenProject';
ref  = fullfile(root, 'rotation_validation_test_results');
umi  = fullfile(root, 'ver_PPP');

modes  = {'simple_shear', 'uniaxial_tension', 'uniaxial_compression', 'plane_strain'};
labels = {'Simple Shear', 'Uniaxial Tension', 'Uniaxial Compression', 'Plane Strain'};
pcts   = [30, 50, 100];

fprintf('\n%-22s %4s %9s %9s %9s\n', 'mode', 'E%', 'TI expl', 'TI impl', 'dTI %');
for i = 1:length(modes)
    for p = pcts
        fe = fullfile(ref, sprintf('orientations_%s_%dpct.txt', modes{i}, p));
        fu = fullfile(umi, sprintf('texture_umat_%s_%dpct.txt', modes{i}, p));
        if ~exist(fe,'file') || ~exist(fu,'file'), continue; end

        de = readmatrix(fe);
        du = readmatrix(fu, 'NumHeaderLines', 1);
        qe = quaternion(de(:,1), de(:,2), de(:,3), de(:,4));
        qu = quaternion(du(:,2), du(:,3), du(:,4), du(:,5));
        ori_e = orientation(qe, cs);
        ori_u = orientation(qu, cs);

        odf_e = calcDensity(ori_e, 'halfwidth', 5*degree);
        odf_u = calcDensity(ori_u, 'halfwidth', 5*degree);
        te = textureindex(odf_e);
        tu = textureindex(odf_u);

        if length(ori_e) == length(ori_u)
            mis = angle(ori_e, ori_u) / degree;
            fprintf('%-22s %4d %9.2f %9.2f %+8.1f%%  miso mean %.1f med %.1f max %.1f deg\n', ...
                labels{i}, p, te, tu, 100*(tu-te)/te, mean(mis), median(mis), max(mis));
        else
            fprintf('%-22s %4d %9.2f %9.2f %+8.1f%%\n', labels{i}, p, te, tu, 100*(tu-te)/te);
        end
    end
end
fprintf('DONE\n');

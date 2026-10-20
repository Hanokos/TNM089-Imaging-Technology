
%% flashNoFlashRun.m
% Flash/no-flash method based on Petschnigg et al. (2004), Section 4.
% Requires Image Processing Toolbox.

close all;
clear;
clc;

imgDir = uigetdir(pwd, 'Select flash/no-flash image folder');
if isequal(imgDir, 0)
    error('No folder selected.');
end

names = {'cave01', 'potsdetail', 'lamp', 'carpet', 'puppets'};

%% Parameters
sdDenoise = 24;
srJoint = 0.001;
srBasic = 0.075;

sdDetail = 12;
srDetail = 0.10;
detailGain = 1;
epsv = 0.02;

tauShad = 0.02;
specLevel = 0.95;
autoExposure = true;
useMask = true;

diffMax = 0.15;

%% Process each image pair
for n = 1:numel(names)

    nm = names{n};

    A = im2double(imread(fullfile(imgDir, ...
        [nm '_01_noflash.tif'])));
    F = im2double(imread(fullfile(imgDir, ...
        [nm '_00_flash.tif'])));

    % Convert grayscale images to RGB if necessary.
    if size(A, 3) == 1
        A = repmat(A, [1 1 3]);
        F = repmat(F, [1 1 3]);
    end

    A = min(max(A, 0), 1);
    F = min(max(F, 0), 1);

    [H, W, ~] = size(A);

    fprintf('\nProcessing %s (%d x %d)...\n', nm, W, H);
    tic;

    % Parameters for the especially noisy cave01 image.
    sdD = sdDenoise;
    sjr = srJoint;
    sdT = sdDetail;
    srT = srDetail;

    if strcmp(nm, 'cave01')
        sdD = 12;
        sjr = 0.05;
        sdT = 12;
        srT = 0.10;
    end

    %% 1. Denoising
    Abase = fastJbf(A, A, sdD, srBasic);
    Anr = fastJbf(A, F, sdD, sjr);

    %% 2. Detail transfer
    Fbase = fastJbf(F, F, sdT, srT);
    Fdetail = ((F + epsv) ./ (Fbase + epsv)) .^ detailGain;

    % Gaussian-based detail layer, used for comparison.
    Fgauss = imgaussfilt(F, sdT);
    DetailG = ((F + epsv) ./ (Fgauss + epsv)) .^ detailGain;

    Atransfer = min(max(Anr .* Fdetail, 0), 1);

    %% 3. Flash-artifact mask
    Alin = A .^ 2.2;
    Flin = F .^ 2.2;

    if autoExposure
        ambientGreen = Alin(:,:,2);
        flashGreen = Flin(:,:,2);

        valid = ambientGreen > 0.02 & flashGreen > 0.02;
        ratios = ambientGreen(valid) ./ flashGreen(valid);

        if ~isempty(ratios)
            ratios = sort(ratios);
            idx = max(1, round(0.95 * numel(ratios)));
            Flin = Flin * ratios(idx);
        end
    end

    shadow = lum(Flin - Alin) <= tauShad;
    shadow = imerode(shadow, strel('disk', 2));
    shadow = imfill(shadow, 'holes');
    shadow = imdilate(shadow, strel('disk', 5));

    spec = lum(F) > specLevel;
    spec = imfill(spec, 'holes');
    spec = imdilate(spec, strel('disk', 5));

    M = imgaussfilt(double(shadow | spec), 3);
    M3 = repmat(M, [1 1 3]);

    fprintf('Mask coverage: %.1f%%\n', 100 * mean(M(:) > 0.5));

    %% Final results
    AfinMask = min(max((1 - M3) .* Anr .* Fdetail ...
        + M3 .* Abase, 0), 1);

    if useMask
        Afin = AfinMask;
    else
        Afin = Atransfer;
    end

    dBJ = mean(abs(Abase - Anr), 3);

    fprintf('Processing completed in %.1f seconds.\n', toc);

    %% Select zoom region
    fsel = figure('Name', 'Select zoom region');
    imshow(Afin);
    title([nm ': Select a detail region and double-click']);

    rz = getrect(gca);
    close(fsel);

    if rz(3) < 3 || rz(4) < 3
        rz = [W*0.4, H*0.4, W*0.15, H*0.15];
    end

    xl = [max(1, rz(1)), min(W, rz(1) + rz(3))];
    yl = [max(1, rz(2)), min(H, rz(2) + rz(4))];

    %% Figure 1: Denoising and detail transfer
    showFig([nm ': Denoising and detail transfer'], ...
        {A, F, Abase, Anr, dBJ, Fdetail/2, Afin}, ...
        {'No-flash image', ...
         'Flash image', ...
         'Basic bilateral filter', ...
         'Flash-guided joint bilateral filter', ...
         sprintf('Basic vs joint difference (mean %.3f)', mean(dBJ(:))), ...
         'Flash detail layer', ...
         'Final result'}, ...
        [2 4], 5, xl, yl, diffMax);

    %% Figure 2: Detail layers and mask comparison
    showFig([nm ': Detail and mask comparison'], ...
        {DetailG/2, Fdetail/2, M, Atransfer, AfinMask, A}, ...
        {'Gaussian-based detail layer', ...
         'Bilateral-based detail layer', ...
         'Flash-artifact mask', ...
         'Detail transfer without mask', ...
         'Detail transfer with mask', ...
         'No-flash reference'}, ...
        [2 3], 0, xl, yl, diffMax);

    drawnow;
end

%% Helper functions
function Y = lum(I)
    Y = 0.299*I(:,:,1) + 0.587*I(:,:,2) + 0.114*I(:,:,3);
end

function showFig(name, imgs, titles, layout, diffIdx, xl, yl, diffMax)

    figure('Name', name, 'NumberTitle', 'off', ...
        'Color', 'w', 'Units', 'normalized', ...
        'Position', [0.02 0.1 0.96 0.75]);

    ax = gobjects(1, numel(imgs));

    for k = 1:numel(imgs)
        ax(k) = subplot(layout(1), layout(2), k);

        if k == diffIdx
            imshow(imgs{k}, [0 diffMax]);
            colormap(ax(k), 'hot');
        else
            imshow(imgs{k});
        end

        title(titles{k}, 'FontSize', 10);
    end

    linkaxes(ax, 'xy');
    xlim(ax(1), xl);
    ylim(ax(1), yl);

    sgtitle(name, 'FontWeight', 'bold');
end

function out = fastJbf(A, G, sd, sr)
    % Fast joint bilateral filter using intensity-level interpolation.

    A = single(A);
    G = single(G);

    sr = max(sr, 0.01);
    nL = min(151, max(6, ceil(2/sr) + 1));
    levels = linspace(0, 1, nL);

    out = zeros(size(A), 'single');

    for k = 1:size(A, 3)
        Ak = A(:,:,k);
        Gk = G(:,:,k);

        pos = Gk * (nL - 1) + 1;
        acc = zeros(size(Ak), 'single');

        for j = 1:nL
            weights = max(0, 1 - abs(pos - j));

            if ~any(weights(:))
                continue;
            end

            rangeWeights = exp(-(Gk - levels(j)).^2 / (2*sr^2));

            numerator = imgaussfilt(rangeWeights .* Ak, sd);
            denominator = imgaussfilt(rangeWeights, sd);

            filtered = numerator ./ max(denominator, 1e-6);
            acc = acc + weights .* filtered;
        end

        out(:,:,k) = acc;
    end

    out = double(out);
end

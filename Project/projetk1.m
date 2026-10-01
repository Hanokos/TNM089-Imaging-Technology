%% flashNoFlashRun.m -- Section 4 of Petschnigg et al. 2004: denoising + detail transfer (+ flash artifact mask)
% Full image resolution, no resizing. Needs Image Processing Toolbox (R2018b+ for sgtitle).
close all; clear; clc;
imgDir = uigetdir(pwd, 'Select the folder with the flash/no-flash images');
if isequal(imgDir, 0), error('No folder selected.'); end

names        = {'lamp','cave01','carpet','puppets','potsdetail'};
pauseBetween = false;   % true: wait for a key press before starting the next scene
closeBetween = false;   % true: close all figures when the next scene starts

%% ---- PARAMETERS (full-resolution pixels) ----
sdDenoise  = 2;      % spatial sigma (paper: 24-48 px neighbourhood on 6MP; smaller here to keep thin detail)
srJoint    = 0.001;  % joint bilateral range sigma, paper value 0.1% (the fast filter floors it at 0.01)
srBasic    = 0.075;  % basic bilateral range sigma, paper 5-10%
sdDetail   = 12;     % spatial sigma of the bilateral filter that makes F_base (larger = more detail transferred)
srDetail   = 0.10;   % range sigma of that filter (larger = more detail, more halos)
detailGain = 1.5;    % 1 = paper; >1 exaggerates the detail layer (our addition)
epsv       = 0.02;   % epsilon in Eq. 6 (paper)
tauShad    = 0.02;   % shadow threshold on linear |flash - ambient|, Eq. 8
specLevel  = 0.95;   % specular = flash luminance above 95% of range (paper)
autoExposure = true; % estimate exposure ratio of the pair (stands in for Eq. 1; needed for the mask)
useMask    = true;   % true: Eq. 5 and 7 with the mask (paper); false: no mask (flash trusted everywhere)
diffMax    = 0.15;   % change that shows as full white in the difference image

for n = 1:numel(names)
    nm = names{n};
    A = im2double(imread(fullfile(imgDir,[nm '_01_noflash.tif'])));   % ambient (no flash)
    F = im2double(imread(fullfile(imgDir,[nm '_00_flash.tif'])));     % flash
    if size(A,3)==1, A = repmat(A,[1 1 3]); F = repmat(F,[1 1 3]); end
    A = min(max(A,0),1);  F = min(max(F,0),1);
    [H,W,~] = size(A);
    fprintf('\n%s: processing full size %d x %d ...\n', nm, W, H);

    sdD = sdDenoise; sjr = srJoint; sdT = sdDetail; srT = srDetail;
    if strcmp(nm,'cave01'), sdD = 12; sjr = 0.05; sdT = 12; srT = 0.10; end  % extremely noisy pair: wide kernels (paper)
    tic;

    %% 4.1 Denoising
    Abase = fastJbf(A, A, sdD, srBasic);     % basic bilateral filter, Eq. 2-3 (also the fallback inside the mask)
    Anr   = fastJbf(A, F, sdD, sjr);         % joint bilateral filter, flash as guide, Eq. 4

    %% 4.2 Flash-to-ambient detail transfer
    Fbase   = fastJbf(F, F, sdT, srT);                         % F_base: basic bilateral filter of the flash image
    Fdetail = ((F+epsv) ./ (Fbase+epsv)) .^ detailGain;        % detail layer (ratio image), Eq. 6
    Fgauss  = imgaussfilt(F, sdT);                             % Gaussian base, only to show the halos (Fig. 5)
    DetailG = ((F+epsv) ./ (Fgauss+epsv)) .^ detailGain;
    Atransfer = min(max(Anr .* Fdetail, 0), 1);                % detail transfer without mask (Fig. 4f)

    %% 4.3 Flash shadows and specularities -> mask M
    Alin = A.^2.2;  Flin = F.^2.2;                             % rough linearization (gamma 2.2)
    if autoExposure   % pixels the flash barely lights (shadows) should satisfy Flin = Alin
        v   = Alin(:,:,2) ./ max(Flin(:,:,2), eps);
        ok0 = (Alin(:,:,2) > 0.02) & (Flin(:,:,2) > 0.02);
        v   = sort(v(ok0));
        Flin = Flin * v(max(1, round(0.95*numel(v))));
    end
    shadow = lum(Flin - Alin) <= tauShad;                                          % Eq. 8
    shadow = imdilate(imfill(imerode(shadow, strel('disk',2)), 'holes'), strel('disk',5));  % erode, fill, dilate
    spec   = imdilate(imfill(lum(F) > specLevel, 'holes'), strel('disk',5));       % saturated in the flash image
    M      = imgaussfilt(double(shadow | spec), 3);            % union, feathered edges
    M3     = repmat(M, [1 1 3]);
    fprintf('  mask covers %.1f%% of the image\n', 100*mean(M(:)>0.5));

    %% Final images (Eq. 5 and 7)
    AnrMask   = (1-M3).*Anr + M3.*Abase;                                          % Eq. 5 (denoising with mask)
    AfinMask  = min(max((1-M3).*Anr.*Fdetail + M3.*Abase, 0), 1);                 % Eq. 7 (with mask)
    Afin      = Atransfer; if useMask, Afin = AfinMask; end                       % our main result
    dBJ       = mean(abs(Abase - Anr), 3);                                        % Fig. 4d: basic vs joint
    fprintf('  done in %.1f s\n', toc);

    %% Pick a zoom region (all panels are linked: zoom one, all follow)
    fsel = figure('Name','Pick'); imshow(Afin);
    title([nm ': DRAG a rectangle for the ZOOM region (something with detail), then double-click']);
    rz = getrect(gca);  close(fsel);
    if isempty(rz) || rz(3)<3, rz = [W*0.4 H*0.4 W*0.15 H*0.15]; end
    xl = [max(1,rz(1)), min(W,rz(1)+rz(3))];  yl = [max(1,rz(2)), min(H,rz(2)+rz(4))];

    %% Figure 1: like the paper's Fig. 4
    showFig([nm ' 1: denoising + detail transfer (paper Fig. 4)'], ...
        {A, F, Abase, Anr, dBJ, Fdetail/2, Afin}, ...
        {'(a) No-flash (ambient)', '(a) Flash', '(b) Denoised: basic bilateral', ...
         '(c) Denoised: joint bilateral (flash-guided)', ...
         sprintf('(d) Difference (b) vs (c), mean %.3f', mean(dBJ(:))), ...
         '(e) Detail layer from flash', '(f) FINAL: denoise + detail transfer'}, ...
        [2 4], 5, xl, yl, diffMax);

    %% Figure 2: bilateral vs Gaussian base (paper Fig. 5) and the mask (paper Fig. 8)
    showFig([nm ' 2: detail layer base (Fig. 5) and flash artifact mask (Fig. 8)'], ...
        {DetailG/2, Fdetail/2, M, Atransfer, AfinMask, A}, ...
        {'Detail layer, Gaussian base (more halos)', 'Detail layer, bilateral base (fewer halos)', ...
         'Mask M (white = flash shadow / specular, flash not trusted)', ...
         'Detail transfer WITHOUT mask', 'Detail transfer WITH mask', 'No-flash (reference)'}, ...
        [2 3], 0, xl, yl, diffMax);
    drawnow;

    if n < numel(names)
        if pauseBetween
            fprintf('Inspect the figures, then press any key in the Command Window for the next scene...\n');
            pause;
        end
        if closeBetween, close all; end
    end
end

%% ============================ helper functions ============================
function Y = lum(I)
Y = 0.299*I(:,:,1) + 0.587*I(:,:,2) + 0.114*I(:,:,3);
end

function showFig(name, imgs, titles, layout, diffIdx, xl, yl, diffMax)
% Shows panels in a grid; panel diffIdx is a difference image (heat map). All axes linked.
figure('Name',name,'NumberTitle','off','Color','w','Units','normalized','Position',[0.02 0.1 0.96 0.75]);
ax = gobjects(1, numel(imgs));
for k = 1:numel(imgs)
    ax(k) = subplot(layout(1), layout(2), k);
    if k == diffIdx, imshow(imgs{k}, [0 diffMax]); colormap(ax(k),'hot');
    else, imshow(imgs{k}); end
    title(titles{k}, 'FontSize', 10);
end
linkaxes(ax, 'xy');  xlim(ax(1), xl);  ylim(ax(1), yl);
sgtitle(name, 'FontWeight', 'bold');
end

function out = fastJbf(A, G, sd, sr)
% Fast joint bilateral filter, piecewise-linear approximation (Durand & Dorsey 2002),
% the speed-up the paper suggests. With G = A this is the basic bilateral filter.
A = single(A);  G = single(G);
sr  = max(sr, 1/100);                      % floor: paper's 0.001 would need ~1000 intensity levels
nL  = min(151, max(6, ceil(2/sr) + 1));    % level spacing <= sr/2
lev = linspace(0, 1, nL);
out = zeros(size(A), 'single');
for k = 1:size(A,3)
    Ak = A(:,:,k);  Gk = G(:,:,k);
    pos = Gk * (nL-1) + 1;
    acc = zeros(size(Ak), 'single');
    for j = 1:nL
        tent = max(0, 1 - abs(pos - j));
        if ~any(tent(:)), continue; end
        Wj  = exp(-(Gk - lev(j)).^2 / (2*sr^2));
        Jj  = imgaussfilt(Wj .* Ak, sd) ./ max(imgaussfilt(Wj, sd), 1e-6);
        acc = acc + tent .* Jj;
    end
    out(:,:,k) = acc;
end
out = double(out);
end
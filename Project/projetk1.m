%% flashNoFlashRun.m -- flash/no-flash, simplified but complete (after Petschnigg et al. 2004)
% Full image resolution, no resizing. Needs Image Processing Toolbox (R2018b+ for sgtitle).
close all; clear; clc;
imgDir = uigetdir(pwd, 'Select the folder with the flash/no-flash images');
if isequal(imgDir, 0), error('No folder selected.'); end

names       = {'lamp','cave01','carpet','puppets','potsdetail'};
stepThrough = true;     % press a key between scenes (figures close)

%% ---- TUNING KNOBS (full-resolution pixels) ----
sdDenoise  = 2;      % spatial sigma; paper ~6-12 gave too much blur on thin strands, so smaller here
srJoint    = 0.001;  % joint range sigma, paper value 0.001 (filter floors it at 0.01, see fastJbf)
srBasic    = 0.075;  % basic bilateral range sigma (paper 5-10%)
sdDetail   = 12;     % spatial sigma, detail layer (larger = more detail transferred)
srDetail   = 0.10;   % range sigma, detail layer   (larger = more detail, more halos)
detailGain = 1.5;    % 1 = paper (no gain); >1 exaggerates detail (our addition, try 1.5-2.5)
epsv       = 0.02;   % epsilon in the detail ratio (paper)
expRatio   = 1;      % (ISO_A*t_A)/(ISO_F*t_F), Eq.1; 1 if unknown
tauShad    = 0.02;   % flash-shadow threshold on linear difference, Eq.8
diffMax    = 0.15;   % change that shows as full white in difference image
useMask    = false;  % false = ignore the flash shadow/specular mask (recommended unless the mask panel looks right)
autoExposure = true;  % estimate the exposure ratio from the images (Eq.1 without camera settings); needed for a sensible mask

for n = 1:numel(names)
    nm = names{n};
    A = im2double(imread(fullfile(imgDir,[nm '_01_noflash.tif'])));   % no-flash (ambient)
    F = im2double(imread(fullfile(imgDir,[nm '_00_flash.tif'])));     % flash
    if size(A,3)==1, A = repmat(A,[1 1 3]); F = repmat(F,[1 1 3]); end
    A = min(max(A,0),1);  F = min(max(F,0),1);
    [H,W,~] = size(A);
    fprintf('\n%s: processing full size %d x %d ...\n', nm, W, H);

    sdD = sdDenoise; sjr = srJoint; sdT = sdDetail; srT = srDetail;
    if strcmp(nm,'cave01'), sdD = 12; sjr = 0.05; sdT = 12; srT = 0.10; end  % very noisy pair
    tic;

    %% Sec. 3: rough linearization (gamma 2.2) + exposure compensation (Eq. 1)
    Alin = A.^2.2;  Flin = (F.^2.2) * expRatio;
    if autoExposure   % pixels the flash barely lights (shadows) should satisfy Flin = Alin
        v   = Alin(:,:,2) ./ max(Flin(:,:,2), eps);
        ok0 = (Alin(:,:,2) > 0.02) & (Flin(:,:,2) > 0.02);
        v   = sort(v(ok0));
        Flin = Flin * v(max(1, round(0.95*numel(v))));
    end

    %% Sec. 4.1: denoising
    Abase = fastJbf(A, A, sdD, srBasic);          % basic bilateral (Eq. 2-3), also the fallback
    Anr   = fastJbf(A, F, sdD, sjr);              % joint bilateral, flash as guide (Eq. 4)

    %% Sec. 4.2: detail transfer
    Fbase   = fastJbf(F, F, sdT, srT);            % bilateral base of flash
    Fdetail = ((F+epsv) ./ (Fbase+epsv)) .^ detailGain;      % detail layer (Eq. 6)
    Fgauss  = imgaussfilt(F, sdT);                % Gaussian base, only to show halos (Fig. 5)
    DetailG = ((F+epsv) ./ (Fgauss+epsv)) .^ detailGain;

    %% Sec. 4.3: flash shadow + specular mask
    shadow = lum(Flin - Alin) <= tauShad;                       % Eq. 8
    shadow = imdilate(imfill(imerode(shadow, strel('disk',2)), 'holes'), strel('disk',5));
    spec   = imdilate(imfill(lum(F) > 0.95, 'holes'), strel('disk',5));   % saturated flash
    M      = imgaussfilt(double(shadow | spec), 3);             % union, feathered
    fprintf('  mask covers %.1f%% of the image (if this is large, the mask is wrong -> keep useMask = false)\n', 100*mean(M(:)>0.5));
    Mshow = M;                        % detected mask, kept for display even if useMask = false
    if ~useMask, M = 0*M; end
    M3     = repmat(M, [1 1 3]);

    %% Final result (Eq. 5 and 7): trust flash only outside the mask
    NoMask = min(max(Anr .* Fdetail, 0), 1);                    % without mask
    Afin   = min(max((1-M3).*Anr.*Fdetail + M3.*Abase, 0), 1);  % with mask = best result
    dif    = mean(abs(Afin - A), 3);
    M3real   = repmat(Mshow, [1 1 3]);                         % the detected mask, used or not
    AfinMask = min(max((1-M3real).*Anr.*Fdetail + M3real.*Abase, 0), 1);   % result WITH mask (Eq. 5, 7)

    %% Sec. 5: white balance from flash (Eq. 9-10)
    D  = Flin - Alin;
    ok = all(Alin >= 0.02, 3) & all(D >= 0.02, 3);
    r  = Alin ./ max(D, eps);  r = r ./ sum(r, 3);              % ambient colour per pixel
    c  = zeros(1,3);  for k = 1:3, t = r(:,:,k); c(k) = mean(t(ok)); end
    c  = c / mean(c);
    Awb = min(A ./ reshape(c,1,1,3), 1);

    %% Sec. 6: continuous flash in YCbCr (Eq. 11), extrapolate Y only
    Cflash = @(a) contFlash(A, F, a);
    fprintf('  done in %.1f s | estimated ambient colour [%.2f %.2f %.2f]\n', toc, c);

    %% Pick regions: one to zoom on, one FLAT area to measure noise
    fsel = figure('Name','Pick'); imshow(Afin);
    title([nm ': DRAG a rectangle for the ZOOM region (something with detail), then double-click']);
    rz = getrect(gca);  close(fsel);
    fsel = figure('Name','Pick'); imshow(Afin);
    title([nm ': DRAG a rectangle on a FLAT area (e.g. plain wall) for the NOISE measurement, then double-click']);
    rf = getrect(gca);  close(fsel);
    if isempty(rz) || rz(3)<3, rz = [W*0.4 H*0.4 W*0.15 H*0.15]; end
    if isempty(rf) || rf(3)<3, rf = [W*0.1 H*0.1 W*0.05 H*0.05]; end
    xl = [max(1,rz(1)), min(W,rz(1)+rz(3))];  yl = [max(1,rz(2)), min(H,rz(2)+rz(4))];

    %% Figures (all panels in a figure are linked: zoom one, all follow)
    showFig([nm ' 1: MAIN comparison'], {A, Afin, dif}, ...
        {'INPUT: No-flash', 'RESULT: Denoised + detail transfer', ...
         sprintf('DIFFERENCE (brighter = bigger change, mean %.3f)', mean(dif(:)))}, ...
        [1 3], 3, xl, yl, diffMax);
    showFig([nm ' 2: how the method works'], ...
        {Abase, Anr, Fdetail/2, DetailG/2, Mshow, NoMask, AfinMask}, ...
        {'Basic bilateral (Fig.4b)', 'Joint bilateral, flash-guided (Fig.4c)', ...
         'Detail layer, bilateral base (Fig.4e)', 'Detail layer, Gaussian base (more halos, Fig.5)', ...
         'Detected flash shadow/specular mask (white = flash not trusted)', 'Result WITHOUT mask', 'Result WITH mask (flash not trusted in white areas)'}, ...
        [2 4], 0, xl, yl, diffMax);
    showFig([nm ' 3: white balance + continuous flash'], ...
        {A, Awb, Cflash(0.5), Cflash(1.5)}, ...
        {'No-flash', 'White balanced (Sec.5)', 'Continuous flash, alpha = 0.5', ...
         'Continuous flash, alpha = 1.5 (extrapolated)'}, [1 4], 0, xl, yl, diffMax);
    drawnow;

    %% Noise + colour numbers (flat patch)
    ix = max(1,round(rf(2))):min(H,round(rf(2)+rf(4)));  jx = max(1,round(rf(1))):min(W,round(rf(1)+rf(3)));
    P = @(I) I(ix, jx, :);
    L = @(I) lum(P(I));
    fprintf('  --- flat patch %d x %d px ---\n', numel(jx), numel(ix));
    fprintf('  noise (std of luminance):  no-flash %.4f | joint bilateral %.4f | final %.4f | flash %.4f\n', ...
            std(L(A),0,'all'), std(L(Anr),0,'all'), std(L(Afin),0,'all'), std(L(F),0,'all'));
    fprintf('  mean colour RGB:  no-flash [%.3f %.3f %.3f] | final [%.3f %.3f %.3f]\n', ...
            squeeze(mean(mean(P(A),1),2)), squeeze(mean(mean(P(Afin),1),2)));

    if stepThrough && n < numel(names)
        fprintf('Inspect the figures, then press any key in the Command Window for the next scene...\n');
        pause; close all;
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

function I = contFlash(A, F, a)
% Eq. 11 in YCbCr; only Y extrapolates, Cb/Cr clamped to the range spanned by A and F.
YA = rgb2ycbcr(A);  YF = rgb2ycbcr(F);
Y  = (1-a)*YA + a*YF;
for k = 2:3
    Y(:,:,k) = min(max(Y(:,:,k), min(YA(:,:,k),YF(:,:,k))), max(YA(:,:,k),YF(:,:,k)));
end
I = min(max(ycbcr2rgb(Y), 0), 1);
end

function out = fastJbf(A, G, sd, sr)
% Fast joint bilateral filter, piecewise-linear approximation (Durand & Dorsey 2002).
% With G = A this is the basic bilateral filter.
A = single(A);  G = single(G);
sr  = max(sr, 1/100);                        % floor: 101 levels max (paper's 0.001 would need ~1000 levels)
nL  = min(151, max(6, ceil(2/sr) + 1));    % spacing <= sr/2 for a smoother approximation
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
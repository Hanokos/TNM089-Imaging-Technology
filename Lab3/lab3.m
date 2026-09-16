%% TNM089 - Imaging Technology
%  Lab 3: High Dynamic Range Imaging
%

clear; close all; clc;

%% ------------------------------------------------------------------
%  Task 1 - The images
%  ------------------------------------------------------------------
pwd
dir('HDR_Images')
imgDir = 'HDR_Images';

% Exposure times, t_j (s). t_{j+1}/t_j = 2
t = [1/4000 1/2000 1/1000 1/500 1/250 1/125 1/60 1/30 1/15 1/6 1/3 2/3 4/3];
N = numel(t);                   % number of exposures (13)

B_row = log2(t);                % log2 exposure time, one value per image

% ---- read the images -------------------------------------------------
% Try Img1.tif / Img1.tiff / img1.tif ...
imgs = cell(1,N);
for j = 1:N
    fname = '';
    cand = {sprintf('Img%d.tif',j),  sprintf('Img%d.tiff',j), ...
            sprintf('img%d.tif',j),  sprintf('img%d.tiff',j), ...
            sprintf('Img%02d.tif',j),sprintf('Img%02d.tiff',j)};
    for c = 1:numel(cand)
        if exist(fullfile(imgDir,cand{c}),'file')
            fname = fullfile(imgDir,cand{c}); break;
        end
    end
    if isempty(fname)
        error('Could not find image %d in folder "%s".', j, imgDir);
    end
    imgs{j} = imread(fname);                 % uint8, 683 x 1024 x 3
end

[rows, cols, ch] = size(imgs{1});
fprintf('Loaded %d images of size %d x %d x %d\n', N, rows, cols, ch);

% ---- display all exposures ------------------------------------------
figure('Name','Task 1: the exposure sequence');
for j = 1:N
    subplot(3,5,j); imshow(imgs{j});
    title(sprintf('Img%d,  t = %.4f s', j, t(j)), 'FontSize', 8);
end

% ---- stack images into one array for convenience ---------------------
% Zstack(i, j, c) = value of pixel i in image j, channel c
Zstack = zeros(rows*cols, N, ch, 'uint8');
for j = 1:N
    im = imgs{j};
    for c = 1:ch
        tmp = im(:,:,c);
        Zstack(:, j, c) = tmp(:);
    end
end

% Intensity (grayscale) version, used for picking the sample pixels
Ystack = zeros(rows*cols, N);
for j = 1:N
    g = rgb2gray(imgs{j});
    Ystack(:,j) = double(g(:));
end

% ---- pick the three pixels -------------------------------------------
% (a) brightest pixel in the darkest image (Img1)
[~, idxBright] = max(Ystack(:,1));
% (b) darkest pixel in the lightest image (Img13)
[~, idxDark]   = min(Ystack(:,N));
% (c) a mid-grey pixel in Img9
[~, idxMid]    = min(abs(Ystack(:,9) - 128));

[rB,cB] = ind2sub([rows cols], idxBright);
[rD,cD] = ind2sub([rows cols], idxDark);
[rM,cM] = ind2sub([rows cols], idxMid);
fprintf('Bright pixel (row,col) = (%d,%d)\n', rB, cB);
fprintf('Dark   pixel (row,col) = (%d,%d)\n', rD, cD);
fprintf('Medium pixel (row,col) = (%d,%d)\n', rM, cM);

% show where they are
figure('Name','Task 1: pixel locations'); imshow(imgs{9}); hold on;
plot(cB,rB,'r+','MarkerSize',12,'LineWidth',2);
plot(cD,rD,'b+','MarkerSize',12,'LineWidth',2);
plot(cM,rM,'g+','MarkerSize',12,'LineWidth',2);
legend('bright','dark','medium'); title('Selected pixels (shown on Img9)');

% ---- the three sequences ---------------------------------------------
seqBright = Ystack(idxBright, :);
seqDark   = Ystack(idxDark,   :);
seqMid    = Ystack(idxMid,    :);

figure('Name','Task 1: pixel value vs exposure');
subplot(1,2,1);
plot(1:N, seqBright,'r-o', 1:N, seqDark,'b-o', 1:N, seqMid,'g-o','LineWidth',1.2);
xlabel('image number j'); ylabel('pixel value z_{ij}');
legend('bright pixel','dark pixel','medium pixel','Location','southeast');
grid on; ylim([0 260]); title('Linear x-axis (image index)');

subplot(1,2,2);
semilogx(t, seqBright,'r-o', t, seqDark,'b-o', t, seqMid,'g-o','LineWidth',1.2);
xlabel('exposure time t_j [s]'); ylabel('pixel value z_{ij}');
grid on; ylim([0 260]); title('Log exposure-time axis');

% Discussion (write this in your report):
%  - Every curve is monotonically increasing: more exposure -> higher value.
%  - The bright pixel saturates (clips at 255) after only a few exposures;
%    all values at 255 are unreliable - the true irradiance is unknown,
%    only a lower bound is known.
%  - The dark pixel stays near 0 for most exposures; those values are
%    dominated by noise/quantisation and are equally unreliable.
%  - Values in the middle of the range (roughly 20-235) are the reliable
%    ones - this is exactly why the weighting function w(z) is a hat
%    function that down-weights values near 0 and 255.
%  - Between saturation and noise the curve is not a straight line, which
%    is the non-linear radiometric response g we estimate in Task 2.


%% ------------------------------------------------------------------
%  Task 2 - Radiometric response function
%  ------------------------------------------------------------------

% ---- weighting function w(z), Eq. 10.7 (hat function) ----------------
Zmin = 0; Zmax = 255; Zmid = 0.5*(Zmin+Zmax);
z = (Zmin:Zmax)';
w = zeros(256,1);
w(z <= Zmid) = z(z <= Zmid) - Zmin;
w(z >  Zmid) = Zmax - z(z >  Zmid);
w = w / max(w);                 % normalise to [0,1] (scale is arbitrary)
w(w < 1e-6) = 1e-6;             % avoid exactly-zero weights in gsolve

figure('Name','Task 2: weighting function');
plot(0:255, w, 'k', 'LineWidth', 1.5); grid on;
xlabel('pixel value z'); ylabel('w(z)'); title('Hat weighting function');

% ---- sample a subset of pixels ---------------------------------------
% Rule of thumb (Debevec): need Npix*(N-1) > 255, and pixels should cover
% the whole intensity range. We sample ~300 pixels spread over intensity.
Npix = 300;
rng(0);                                        % reproducible sampling
[~, order] = sort(Ystack(:,7));                % sort by mid-exposure value
sampleIdx  = order(round(linspace(1, numel(order), Npix)));

% ---- build Z and B and solve for each channel ------------------------
lambda = 20;                                   % smoothness. Try 1 ... 100

B = repmat(B_row, Npix, 1);   % <-- matrix, because gsolve uses B(i,j)

gfun = zeros(256, ch);        % g for R, G, B
lE   = zeros(Npix, ch);       % log irradiance of the sampled pixels

for c = 1:ch
    Z = double(Zstack(sampleIdx, :, c));       % Npix x N
    [gc, lEc] = gsolve(Z, B, lambda, w);
    gfun(:,c) = gc;
    lE(:,c)   = lEc;
end

% ---- plot g(z) -------------------------------------------------------
figure('Name','Task 2: response function');
subplot(1,2,1);
plot(gfun(:,1), 0:255, 'r', gfun(:,2), 0:255, 'g', gfun(:,3), 0:255, 'b', ...
     'LineWidth', 1.5);
xlabel('log_2 exposure,  g(z)'); ylabel('pixel value z');
legend('R','G','B','Location','southeast'); grid on;
title(sprintf('g(z),  \\lambda = %g', lambda));

% ---- plot 2.^g : mapping pixel value -> irradiance --------------------
subplot(1,2,2);
plot(0:255, 2.^gfun(:,1), 'r', 0:255, 2.^gfun(:,2), 'g', ...
     0:255, 2.^gfun(:,3), 'b', 'LineWidth', 1.5);
xlabel('pixel value z'); ylabel('relative irradiance E = 2^{g(z)}');
legend('R','G','B','Location','northwest'); grid on;
title('Exposure (inverse response) curve');

% Tip: also try set(gca,'YScale','log') on the right plot - on a log axis
% the curve is simply g(z) again, which makes the shape easier to read.


%% ------------------------------------------------------------------
%  Task 3 - Constructing the radiance map
%  ------------------------------------------------------------------

HDR = makeRadianceMap(Zstack, B_row, gfun, w, rows, cols);

% ---- dynamic range ---------------------------------------------------
Yhdr = 0.2126*HDR(:,:,1) + 0.7152*HDR(:,:,2) + 0.0722*HDR(:,:,3);
vals = Yhdr(isfinite(Yhdr) & Yhdr > 0);
DRratio  = max(vals) / min(vals);
DRstops  = log2(DRratio);
DRdB     = 20*log10(DRratio);
fprintf('\nHDR dynamic range: %.3e : 1   (%.1f stops / %.1f dB)\n', ...
        DRratio, DRstops, DRdB);

% A more robust number, ignoring outlier pixels (1st / 99th percentile):
p = prctile(vals, [1 99]);
fprintf('Robust (1%%-99%%) range: %.3e : 1  (%.1f stops)\n', ...
        p(2)/p(1), log2(p(2)/p(1)));
% For comparison, a single 8-bit exposure has about 255:1 = 8 stops.

% ---- pseudo-colour display of the irradiance -------------------------
figure('Name','Task 3: radiance map');
subplot(1,2,1);
imagesc(Yhdr); axis image off; colormap(gca, parula); colorbar;
title('Irradiance, linear scale');
subplot(1,2,2);
imagesc(log2(Yhdr)); axis image off; colormap(gca, parula); colorbar;
title('log_2 irradiance');

% Save the result if you like:
% hdrwrite(HDR, 'result.hdr');


%% ------------------------------------------------------------------
%  Task 4 - Tone mapping
%  ------------------------------------------------------------------

HDRn = HDR / max(HDR(:));       % normalised to [0,1]

% ---- (a) linear scaling ----------------------------------------------
figure('Name','Task 4: linear scaling');
ks = [1 5 20 100];
for n = 1:4
    subplot(2,2,n); imshow(HDRn * ks(n));
    title(sprintf('linear, k = %d', ks(n)));
end
% Why it looks bad: a few very bright pixels (light sources, specular
% highlights) set max(HDR(:)), so everything else is squeezed into the
% lowest few code values -> almost black. Raising k brightens the shadows
% but immediately clips the highlights: a linear map cannot show a
% ~1e5:1 range on an 8-bit, ~100:1 display.

% ---- (b) gamma mapping -----------------------------------------------
figure('Name','Task 4: gamma mapping');
gammas = [1/1.5 1/2.2 1/3 1/5];
for n = 1:4
    subplot(2,2,n); imshow(HDRn .^ gammas(n));
    title(sprintf('gamma = 1/%.1f', 1/gammas(n)));
end

% ---- (c) logarithmic mapping -----------------------------------------
eps0 = 1e-6;
LOG  = log(1 + HDRn/eps0);
LOG  = LOG / max(LOG(:));
figure('Name','Task 4: log mapping'); imshow(LOG);
title('Global log transform');

% Both gamma and log are *global* operators: one curve for all pixels.
% They compress the range, but colours desaturate and local contrast is
% lost, because the same strong compression is applied everywhere.

% ---- (d) MATLAB tonemap ----------------------------------------------
figure('Name','Task 4: tonemap');
subplot(2,2,1); imshow(tonemap(HDR));
title('tonemap, defaults');
subplot(2,2,2); imshow(tonemap(HDR,'AdjustLightness',[0.1 1],'AdjustSaturation',1.5));
title('lighter, more saturated');
subplot(2,2,3); imshow(tonemap(HDR,'AdjustLightness',[0 0.9],'AdjustSaturation',0.8));
title('darker, less saturated');
subplot(2,2,4); imshow(tonemap(HDR,'NumberOfTiles',[8 8]));
title('8x8 tiles');

% ---- (e) localtonemap (local operator, needs single) ------------------
HDRs = single(HDR);
HDRs(HDRs <= 0) = eps('single');
figure('Name','Task 4: localtonemap');
subplot(1,3,1); imshow(localtonemap(HDRs));
title('localtonemap, defaults');
subplot(1,3,2); imshow(localtonemap(HDRs,'RangeCompression',0.3,'EnhanceContrast',0.3));
title('mild compression');
subplot(1,3,3); imshow(localtonemap(HDRs,'RangeCompression',0.9,'EnhanceContrast',0.7));
title('strong compression + contrast');

% localtonemap separates a base layer (compressed) from a detail layer
% (preserved), so local contrast and texture survive - that is why it
% usually looks noticeably better than the global operators above.


%% ------------------------------------------------------------------
%  Local function: radiance map (Eq. 10.11)
%  ------------------------------------------------------------------
function HDR = makeRadianceMap(Zstack, B_row, gfun, w, rows, cols)
% Zstack : (rows*cols) x Nimages x 3 uint8 pixel values
% B_row  : 1 x Nimages, log2 of the exposure times
% gfun   : 256 x 3, log2 exposure response for R,G,B
% w      : 256 x 1 weighting function
%
%            sum_j w(z_ij) * ( g(z_ij) - log t_j )
% log E_i = ---------------------------------------      (Eq. 10.11)
%                    sum_j w(z_ij)

    [Npx, N, ch] = size(Zstack);
    HDR = zeros(rows, cols, ch);

    for c = 1:ch
        num = zeros(Npx,1);
        den = zeros(Npx,1);

        Zc = double(Zstack(:,:,c));
        gc = gfun(:,c);

        for j = 1:N
            idx  = Zc(:,j) + 1;              % 1..256
            wij  = w(idx);
            num  = num + wij .* (gc(idx) - B_row(j));
            den  = den + wij;
        end

        logE = num ./ den;

        % Fallback for pixels that are saturated or black in every image
        % (den ~ 0): use the single exposure whose value is closest to the
        % middle of the range, without weighting.
        bad = (den < 1e-5) | ~isfinite(logE);
        if any(bad)
            [~, jbest] = min(abs(Zc(bad,:) - 128), [], 2);
            zb  = Zc(bad,:);
            lin = sub2ind(size(zb), (1:size(zb,1))', jbest);
            logE(bad) = gc(zb(lin) + 1) - B_row(jbest)';
        end

        HDR(:,:,c) = reshape(2.^logE, rows, cols);   % base-2, since g is log2
    end
end
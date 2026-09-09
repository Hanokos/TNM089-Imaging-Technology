%% Load data
clear all
close all
clc
load('Lab2_data.mat')  % loads MS_Image, Camera, wl

[rows, cols, bands] = size(MS_Image);  % 1040 x 1392 x 31

%% (Optional) plot camera sensitivity functions
figure;
plot(wl, Camera(1,:), 'r', wl, Camera(2,:), 'g', wl, Camera(3,:), 'b', 'LineWidth', 1.5)
xlabel('Wavelength (nm)')
ylabel('Sensitivity')
title('Camera RGB Sensitivity Functions')
legend('R','G','B')
grid on

%% Compute raw sensor RGB image
% Reshape MS_Image so each pixel's spectrum is a row: (rows*cols) x bands
I_reshaped = reshape(MS_Image, rows*cols, bands);   % (1040*1392) x 31

% Multiply by Camera' (31 x 3) -> each column n of Camera' gives
% an = sum(Mn .* I) for every pixel simultaneously
RGB_reshaped = I_reshaped * Camera';                 % (1040*1392) x 3

% Reshape back to image form
RGB_image = reshape(RGB_reshaped, rows, cols, 3);    % 1040 x 1392 x 3

%% Display the raw (linear) sensor image
figure;
imshow(RGB_image)
title('Raw linear sensor RGB image')

%% ========================================================================
%% STEP 1: Create the Bayer-pattern CFA image
%% ========================================================================
% Standard Bayer RGGB pattern:
%   R G R G ...
%   G B G B ...
% (row 1,3,5... : R G R G ; row 2,4,6... : G B G B)

R = RGB_image(:,:,1);
G = RGB_image(:,:,2);
B = RGB_image(:,:,3);

% Create binary masks for each channel
maskR = zeros(rows, cols);
maskG = zeros(rows, cols);
maskB = zeros(rows, cols);

maskR(1:2:end, 1:2:end) = 1;     % R at odd rows, odd cols
maskB(2:2:end, 2:2:end) = 1;     % B at even rows, even cols
maskG = 1 - maskR - maskB;       % G everywhere else (50%)

% CFA image: single-channel mosaic, and a "sparse" RGB representation
CFA = R.*maskR + G.*maskG + B.*maskB;

CFA_RGB = zeros(rows, cols, 3);   % useful for visualizing/processing
CFA_RGB(:,:,1) = R.*maskR;
CFA_RGB(:,:,2) = G.*maskG;
CFA_RGB(:,:,3) = B.*maskB;

figure; imshow(CFA); title('Raw CFA (mosaic) image')

%% ========================================================================
%% METHOD 1: Bilinear interpolation, each channel independently
%% ========================================================================
% R and B are sampled on a grid with period 2 in both directions
% -> use kernel h1
h1 = [1 2 1; 2 4 2; 1 2 1]/4;

% G is sampled in a quincunx (checkerboard) pattern
% -> use kernel h2
h2 = [0 1 0; 1 4 1; 0 1 0]/4;

R1 = conv2(CFA_RGB(:,:,1), h1, 'same');
G1 = conv2(CFA_RGB(:,:,2), h2, 'same');
B1 = conv2(CFA_RGB(:,:,1)*0 + CFA_RGB(:,:,3), h1, 'same'); % B channel

RGB_bilinear = cat(3, R1, G1, B1);
RGB_bilinear = max(min(RGB_bilinear,1),0); % clip to valid range

figure; imshow(RGB_bilinear); title('Method 1: Independent bilinear interpolation')

%% ========================================================================
%% METHOD 2: Bilinear G + spectral correlation (color-difference) for R,B
%% ========================================================================
% Idea: interpolate G fully first (assume it's the "sharpest"/most sampled
% channel). Then compute R-G and B-G only where R (resp. B) is known,
% interpolate these DIFFERENCE images (smoother, low frequency), and add
% G back. This exploits high correlation between color channels.

G2 = G1; % reuse full bilinear G from Method 1

% --- Red ---
diffR = zeros(rows, cols);
diffR(maskR==1) = R(maskR==1) - G2(maskR==1);
diffR_interp = conv2(diffR, h1, 'same');   % same sampling grid as R
R2 = G2 + diffR_interp;

% --- Blue ---
diffB = zeros(rows, cols);
diffB(maskB==1) = B(maskB==1) - G2(maskB==1);
diffB_interp = conv2(diffB, h1, 'same');
B2 = G2 + diffB_interp;

RGB_colordiff = cat(3, R2, G2, B2);
RGB_colordiff = max(min(RGB_colordiff,1),0);

figure; imshow(RGB_colordiff); title('Method 2: Bilinear G + color-difference R,B')

%% ========================================================================
%% METHOD 3: Edge-directed interpolation for G + color-difference for R,B
%% ========================================================================
% Classic Hamilton-Adams style edge-directed interpolation for G:
% At each missing-G location, estimate horizontal & vertical gradients
% using both G and the co-located R or B channel (2nd derivative term),
% then interpolate along the direction of lower gradient.

Gfull = CFA_RGB(:,:,2);   % sparse G (zeros where missing)
G3 = Gfull;

% Pad arrays to avoid boundary issues
pad = 2;
Rp = padarray(R.*maskR + B.*maskB.*0, [pad pad], 'replicate'); % placeholder, replaced below
Gp = padarray(Gfull, [pad pad], 'replicate');
CFAp = padarray(CFA, [pad pad], 'replicate');
maskRp = padarray(maskR, [pad pad], 'replicate');
maskBp = padarray(maskB, [pad pad], 'replicate');

for i = 1+pad:rows+pad
    for j = 1+pad:cols+pad
        if maskRp(i,j)==0 && maskBp(i,j)==0
            continue % this is already a G pixel, skip
        end
        
        % Use the CFA value at this pixel (R or B) for 2nd derivative term
        Cval = CFAp(i,j);
        
        % Horizontal gradient estimate
        gradH = abs(Gp(i,j-1) - Gp(i,j+1)) + ...
                abs(2*Cval - CFAp(i,j-2) - CFAp(i,j+2));
        % Vertical gradient estimate
        gradV = abs(Gp(i-1,j) - Gp(i+1,j)) + ...
                abs(2*Cval - CFAp(i-2,j) - CFAp(i+2,j));
        
        Hinterp = (Gp(i,j-1) + Gp(i,j+1))/2 + (2*Cval - CFAp(i,j-2) - CFAp(i,j+2))/4;
        Vinterp = (Gp(i-1,j) + Gp(i+1,j))/2 + (2*Cval - CFAp(i-2,j) - CFAp(i+2,j))/4;
        
        if gradH < gradV
            Gval = Hinterp;
        elseif gradV < gradH
            Gval = Vinterp;
        else
            Gval = (Hinterp + Vinterp)/2;
        end
        
        G3(i-pad, j-pad) = Gval;
    end
end

G3 = max(min(G3,1),0);

% --- Now use color-difference method for R and B, based on edge-directed G ---
diffR3 = zeros(rows, cols);
diffR3(maskR==1) = R(maskR==1) - G3(maskR==1);
diffR3_interp = conv2(diffR3, h1, 'same');
R3 = G3 + diffR3_interp;

diffB3 = zeros(rows, cols);
diffB3(maskB==1) = B(maskB==1) - G3(maskB==1);
diffB3_interp = conv2(diffB3, h1, 'same');
B3 = G3 + diffB3_interp;

RGB_edgedir = cat(3, R3, G3, B3);
RGB_edgedir = max(min(RGB_edgedir,1),0);

figure; imshow(RGB_edgedir); title('Method 3: Edge-directed G + color-difference R,B')

%% ========================================================================
%% EVALUATION: MSE for each method
%% ========================================================================
MSE_bilinear  = mean((RGB_bilinear(:)  - RGB_image(:)).^2);
MSE_colordiff = mean((RGB_colordiff(:) - RGB_image(:)).^2);
MSE_edgedir   = mean((RGB_edgedir(:)   - RGB_image(:)).^2);

fprintf('MSE (independent bilinear):        %.6e\n', MSE_bilinear);
fprintf('MSE (color-difference):            %.6e\n', MSE_colordiff);
fprintf('MSE (edge-directed + colordiff):   %.6e\n', MSE_edgedir);

%% Visualize difference (error) images to find problem regions
err_bilinear  = sqrt(sum((RGB_bilinear  - RGB_image).^2, 3));
err_colordiff = sqrt(sum((RGB_colordiff - RGB_image).^2, 3));
err_edgedir   = sqrt(sum((RGB_edgedir   - RGB_image).^2, 3));

figure; imshow(RGB_image)
title('Click two corners to define a crop region, then press Enter')
[x, y] = ginput(2);
r = round(min(y)):round(max(y));
c = round(min(x)):round(max(x));

figure;
subplot(2,2,1); imshow(RGB_image(r,c,:)); title('Reference')
subplot(2,2,2); imshow(RGB_bilinear(r,c,:)); title('Bilinear')
subplot(2,2,3); imshow(RGB_colordiff(r,c,:)); title('Color-diff')
subplot(2,2,4); imshow(RGB_edgedir(r,c,:)); title('Edge-directed')

%% ========================================================================
%% TASK 3: Automatic White Balance
%% ========================================================================
% Use the demosaiced RGB image (edge-directed, best quality) as input
RGB_in = RGB_edgedir;   % you can also try this with RGB_bilinear/RGB_colordiff

R = RGB_in(:,:,1);
G = RGB_in(:,:,2);
B = RGB_in(:,:,3);

%% ------------------------------------------------------------------
%% METHOD 1: Gray World Assumption
%% ------------------------------------------------------------------
% Assumption: the average reflectance in a scene, averaged over all
% wavelengths, is achromatic (gray). So the average R, G, B should
% be equal. We compute per-channel gain factors to enforce this,
% normalized so G stays (roughly) unchanged as a reference channel.

meanR = mean(R(:));
meanG = mean(G(:));
meanB = mean(B(:));

% Scale factors: bring each channel's mean to match the G-channel mean
kR_gw = meanG / meanR;
kG_gw = 1;               % reference channel, unchanged
kB_gw = meanG / meanB;

RGB_grayworld = zeros(size(RGB_in));
RGB_grayworld(:,:,1) = R * kR_gw;
RGB_grayworld(:,:,2) = G * kG_gw;
RGB_grayworld(:,:,3) = B * kB_gw;

RGB_grayworld = max(min(RGB_grayworld,1),0);  % clip to valid range

fprintf('Gray World gains: R=%.4f, G=%.4f, B=%.4f\n', kR_gw, kG_gw, kB_gw);

%% ------------------------------------------------------------------
%% METHOD 2: White Patch (Max-RGB) Assumption
%% ------------------------------------------------------------------
% Assumption: the brightest pixel/highlight in the scene is a specular
% reflection of the light source, and should therefore be white/neutral.
% We normalize each channel by its own maximum value.

maxR = max(R(:));
maxG = max(G(:));
maxB = max(B(:));

kR_wp = maxG / maxR;
kG_wp = 1;
kB_wp = maxG / maxB;

RGB_whitepatch = zeros(size(RGB_in));
RGB_whitepatch(:,:,1) = R * kR_wp;
RGB_whitepatch(:,:,2) = G * kG_wp;
RGB_whitepatch(:,:,3) = B * kB_wp;

RGB_whitepatch = max(min(RGB_whitepatch,1),0);

fprintf('White Patch gains: R=%.4f, G=%.4f, B=%.4f\n', kR_wp, kG_wp, kB_wp);

%% ------------------------------------------------------------------
%% METHOD 3: Physical White Reference (Barium Sulfate patch)
%% ------------------------------------------------------------------
% Use the known white reference sub-image to estimate the light source
% color directly, since (by definition) it should reflect all
% wavelengths equally -> its RGB response directly tells us the
% color of the illuminant.

whiteRef = RGB_in(750:850, 1200:1300, :);

refR = mean(mean(whiteRef(:,:,1)));
refG = mean(mean(whiteRef(:,:,2)));
refB = mean(mean(whiteRef(:,:,3)));

kR_ref = refG / refR;
kG_ref = 1;
kB_ref = refG / refB;

RGB_refwhite = zeros(size(RGB_in));
RGB_refwhite(:,:,1) = R * kR_ref;
RGB_refwhite(:,:,2) = G * kG_ref;
RGB_refwhite(:,:,3) = B * kB_ref;

RGB_refwhite = max(min(RGB_refwhite,1),0);

fprintf('Reference White gains: R=%.4f, G=%.4f, B=%.4f\n', kR_ref, kG_ref, kB_ref);

%% ------------------------------------------------------------------
%% DISPLAY COMPARISON
%% ------------------------------------------------------------------
figure;
subplot(2,2,1); imshow(RGB_in);          title('Original (no AWB)')
subplot(2,2,2); imshow(RGB_grayworld);   title('Gray World AWB')
subplot(2,2,3); imshow(RGB_whitepatch);  title('White Patch AWB')
subplot(2,2,4); imshow(RGB_refwhite);    title('Physical White Reference AWB')

%% Highlight where the white reference patch is located (for verification)
figure; imshow(RGB_in); hold on
rectangle('Position', [1200, 750, 100, 100], 'EdgeColor', 'y', 'LineWidth', 2)
title('Location of white reference patch used')

%% Verify: check the corrected color of the white patch itself under each method
patch_gw  = RGB_grayworld(750:850, 1200:1300, :);
patch_wp  = RGB_whitepatch(750:850, 1200:1300, :);
patch_ref = RGB_refwhite(750:850, 1200:1300, :);

fprintf('\nMean RGB of white patch AFTER correction:\n');
fprintf('Gray World:  R=%.4f G=%.4f B=%.4f\n', mean(mean(patch_gw(:,:,1))), mean(mean(patch_gw(:,:,2))), mean(mean(patch_gw(:,:,3))));
fprintf('White Patch: R=%.4f G=%.4f B=%.4f\n', mean(mean(patch_wp(:,:,1))), mean(mean(patch_wp(:,:,2))), mean(mean(patch_wp(:,:,3))));
fprintf('Ref White:   R=%.4f G=%.4f B=%.4f\n', mean(mean(patch_ref(:,:,1))), mean(mean(patch_ref(:,:,2))), mean(mean(patch_ref(:,:,3))));

%% ========================================================================
%% TASK 4: Image Space Rendering — each step shown full-size
%% ========================================================================
RGB_wb = RGB_refwhite;   % starting point from Task 3

% Quick sanity check on the input range
fprintf('Input (RGB_wb): min=%.4f mean=%.4f max=%.4f\n', ...
    min(RGB_wb(:)), mean(RGB_wb(:)), max(RGB_wb(:)))

figure; imshow(RGB_wb); title('Step 0: White-balanced (linear, no processing)')

%% ------------------------------------------------------------------
%% STEP 1: Contrast stretching
%% ------------------------------------------------------------------
low_pct  = 1;
high_pct = 99;

vals = RGB_wb(:);
lo = prctile(vals, low_pct);
hi = prctile(vals, high_pct);

RGB_stretch = (RGB_wb - lo) / (hi - lo);
RGB_stretch = max(min(RGB_stretch, 1), 0);

fprintf('After stretch: min=%.4f mean=%.4f max=%.4f\n', ...
    min(RGB_stretch(:)), mean(RGB_stretch(:)), max(RGB_stretch(:)))

figure; imshow(RGB_stretch); title('Step 1: Contrast stretch')

%% ------------------------------------------------------------------
%% STEP 2: Gamma correction
%% ------------------------------------------------------------------
gamma = 1.8;   % adjust based on how Step 1 looks -- try 1.5 / 1.8 / 2.4 and compare
RGB_gamma = RGB_stretch .^ (1/gamma);

fprintf('After gamma (%.1f): min=%.4f mean=%.4f max=%.4f\n', ...
    gamma, min(RGB_gamma(:)), mean(RGB_gamma(:)), max(RGB_gamma(:)))

figure; imshow(RGB_gamma); title(sprintf('Step 2: Gamma correction (\\gamma=%.1f)', gamma))

%% ------------------------------------------------------------------
%% STEP 3: S-curve tone mapping
%% ------------------------------------------------------------------
intensity = 6;

RGB_s = 1 ./ (1 + exp(-intensity * (RGB_gamma - 0.5)));
s_lo = 1 / (1 + exp(-intensity * (0 - 0.5)));
s_hi = 1 / (1 + exp(-intensity * (1 - 0.5)));
RGB_s = (RGB_s - s_lo) / (s_hi - s_lo);
RGB_s = max(min(RGB_s, 1), 0);

figure; imshow(RGB_s); title(sprintf('Step 3: S-curve tone map (intensity=%d)', intensity))

%% ------------------------------------------------------------------
%% STEP 4: Unsharp masking
%% ------------------------------------------------------------------
sigma_blur = 2;
amount     = 0.8;

RGB_blur = imgaussfilt(RGB_s, sigma_blur);
RGB_highpass = RGB_s - RGB_blur;
RGB_sharp = RGB_s + amount * RGB_highpass;
RGB_sharp = max(min(RGB_sharp, 1), 0);

figure; imshow(RGB_sharp); title(sprintf('Step 4: Unsharp masking (\\sigma=%d, amount=%.1f)', sigma_blur, amount))

%% ------------------------------------------------------------------
%% STEP 5: Convert to uint8 (final output)
%% ------------------------------------------------------------------
RGB_final_uint8 = im2uint8(RGB_sharp);

figure; imshow(RGB_final_uint8); title('Step 5: Final image (uint8)')

%% ------------------------------------------------------------------
%% Final comparison: Task 1 raw reference vs Task 4 final
%% ------------------------------------------------------------------
figure; imshow(RGB_image); title('Task 1: Raw linear reference (for comparison)')
figure; imshow(RGB_final_uint8); title('Task 4: Final rendered image (for comparison)')
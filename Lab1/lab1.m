clear;
close all;
clc;

%% TNM089 - Imaging Technology
% Lab 1 - Autofocus
%
% Task 1:
%   - Select a focus window using ginput
%   - Evaluate VAR, EIG and FT2
%   - Find the best-focused image
%
% Task 2:
%   - Evaluate the focus measures for different window sizes and positions
%   - Plot performance curves
%   - Evaluate unimodality and monotonicity
%   - Compare computational efficiency
%
% Task 3:
%   - Construct a full-focus image using VAR
%   - Compare different focus-window sizes in image quality and time
%
% The implementation below is intended for the KaktusFocus image stack.

%% Paths and image loading

% Keep this if TNM089 is needed for functions/files in your setup.
if isfolder('C:\TNM089')
    addpath(genpath('C:\TNM089'));
end

imageFolder = 'KaktusFocus';

files = dir(fullfile(imageFolder, '*.jpg'));

if isempty(files)
    error(['No JPG images were found in "%s". Make sure the folder ', ...
           'is in the current MATLAB folder or change imageFolder.'], ...
           imageFolder);
end

% Sort files by name so that image numbering is reproducible.
[~, order] = sort({files.name});
files = files(order);

images = cell(1, length(files));

for i = 1:length(files)

    I = imread(fullfile(imageFolder, files(i).name));

    % KaktusFocus should be grayscale, but convert RGB images if needed.
    if ndims(I) == 3
        I = rgb2gray(I);
    end

    images{i} = double(I);

end

N = length(images);

fprintf('Number of images loaded: %d\n', N);

if N ~= 26
    warning(['The lab description states that KaktusFocus contains 26 ', ...
             'grayscale images. The program found %d images.'], N);
end

[height, width] = size(images{1});

fprintf('Image size: %d x %d pixels\n\n', height, width);

%% ================================================================
%  TASK 1 - AUTOFOCUS
%  ================================================================

% Display one image and let the user select two opposite corners
% of the focus window.
figure;
imshow(images{13}, []);
title('Task 1 - Click two opposite corners of the focus window');

[x, y] = ginput(2);

x1 = round(min(x));
x2 = round(max(x));
y1 = round(min(y));
y2 = round(max(y));

% Make sure the selected coordinates are inside the image.
x1 = max(1, min(width, x1));
x2 = max(1, min(width, x2));
y1 = max(1, min(height, y1));
y2 = max(1, min(height, y2));

if x1 == x2 || y1 == y2
    error('The selected focus window is too small. Please run the program again.');
end

fprintf('Selected focus window: x = %d:%d, y = %d:%d\n', ...
        x1, x2, y1, y2);

% Allocate focus-measure vectors.
VAR = zeros(1, N);
EIG = zeros(1, N);
FT2 = zeros(1, N);

% Calculate the three focus measures for every image.
for i = 1:N

    I = images{i}(y1:y2, x1:x2);

    % ------------------------------------------------------------
    % VAR - statistical focus measure
    % ------------------------------------------------------------
    VAR(i) = var(I(:));

    % ------------------------------------------------------------
    % EIG - spatial/filter-based focus measure
    % ------------------------------------------------------------
    Gx = imfilter(I, [-1 0 1], 'replicate');
    Gy = imfilter(I, [-1; 0; 1], 'replicate');

    EIG(i) = sum(Gx(:).^2 + Gy(:).^2);

    % ------------------------------------------------------------
    % FT2 - Fourier-domain focus measure
    % ------------------------------------------------------------
    F = fft2(I);
    FT2(i) = sum(abs(F(:)));

end

% Find the best-focused image according to each measure.
[~, bestVAR] = max(VAR);
[~, bestEIG] = max(EIG);
[~, bestFT2] = max(FT2);

fprintf('\nTask 1 - Best-focused images\n');
fprintf('VAR = image %d (%s)\n', bestVAR, files(bestVAR).name);
fprintf('EIG = image %d (%s)\n', bestEIG, files(bestEIG).name);
fprintf('FT2 = image %d (%s)\n', bestFT2, files(bestFT2).name);

% Plot the focus curves.
figure;

plot(1:N, VAR ./ max(VAR), 'o-', 'LineWidth', 1.2);
hold on;
plot(1:N, EIG ./ max(EIG), 's-', 'LineWidth', 1.2);
plot(1:N, FT2 ./ max(FT2), '^-', 'LineWidth', 1.2);

xlabel('Image number');
ylabel('Normalized focus measure');
title('Task 1 - Focus measures');
legend('VAR', 'EIG', 'FT2', 'Location', 'best');
grid on;

% Show the best image found by each method.
figure;

subplot(1, 3, 1);
imshow(images{bestVAR}, []);
title(sprintf('VAR - image %d', bestVAR));

subplot(1, 3, 2);
imshow(images{bestEIG}, []);
title(sprintf('EIG - image %d', bestEIG));

subplot(1, 3, 3);
imshow(images{bestFT2}, []);
title(sprintf('FT2 - image %d', bestFT2));

sgtitle('Task 1 - Best-focused images');

%% ================================================================
%  TASK 2 - EVALUATION OF FOCUS MEASURES
%  ================================================================

% Focus-window sizes.
windowSizes = [32 64 128];

% Positions are defined as normalized locations in the image:
%   1 = upper-left
%   2 = center
%   3 = lower-right
%
% The positions are calculated for every window size so that the
% entire window remains inside the image. This makes the comparison
% between positions fair.
positionNames = {'Upper-left', 'Center', 'Lower-right'};
positionFractions = [0.25 0.25;
                     0.50 0.50;
                     0.75 0.75];

numSizes = length(windowSizes);
numPositions = size(positionFractions, 1);

resultsVAR = cell(numSizes, numPositions);
resultsEIG = cell(numSizes, numPositions);
resultsFT2 = cell(numSizes, numPositions);

timeVAR = zeros(numSizes, numPositions);
timeEIG = zeros(numSizes, numPositions);
timeFT2 = zeros(numSizes, numPositions);

bestImagesVAR = zeros(numSizes, numPositions);
bestImagesEIG = zeros(numSizes, numPositions);
bestImagesFT2 = zeros(numSizes, numPositions);

% Unimodality and monotonicity results.
unimodalVAR = false(numSizes, numPositions);
unimodalEIG = false(numSizes, numPositions);
unimodalFT2 = false(numSizes, numPositions);

monoVAR = zeros(numSizes, numPositions);
monoEIG = zeros(numSizes, numPositions);
monoFT2 = zeros(numSizes, numPositions);

for w = 1:numSizes

    windowSize = windowSizes(w);

    for r = 1:numPositions

        % Calculate a valid center position for this window size.
        centerX = round(positionFractions(r, 1) * width);
        centerY = round(positionFractions(r, 2) * height);

        x1 = round(centerX - windowSize/2);
        y1 = round(centerY - windowSize/2);

        % Force the complete focus window to stay inside the image.
        x1 = max(1, min(width - windowSize + 1, x1));
        y1 = max(1, min(height - windowSize + 1, y1));

        x2 = x1 + windowSize - 1;
        y2 = y1 + windowSize - 1;

        VAR_temp = zeros(1, N);
        EIG_temp = zeros(1, N);
        FT2_temp = zeros(1, N);

        % --------------------------------------------------------
        % VAR
        % --------------------------------------------------------
        tic;

        for i = 1:N

            I = images{i}(y1:y2, x1:x2);
            VAR_temp(i) = var(I(:));

        end

        timeVAR(w, r) = toc;

        % --------------------------------------------------------
        % EIG
        % --------------------------------------------------------
        tic;

        for i = 1:N

            I = images{i}(y1:y2, x1:x2);

            Gx = imfilter(I, [-1 0 1], 'replicate');
            Gy = imfilter(I, [-1; 0; 1], 'replicate');

            EIG_temp(i) = sum(Gx(:).^2 + Gy(:).^2);

        end

        timeEIG(w, r) = toc;

        % --------------------------------------------------------
        % FT2
        % --------------------------------------------------------
        tic;

        for i = 1:N

            I = images{i}(y1:y2, x1:x2);

            F = fft2(I);
            FT2_temp(i) = sum(abs(F(:)));

        end

        timeFT2(w, r) = toc;

        % Store results.
        resultsVAR{w, r} = VAR_temp;
        resultsEIG{w, r} = EIG_temp;
        resultsFT2{w, r} = FT2_temp;

        % Best-focused images.
        [~, bestImagesVAR(w, r)] = max(VAR_temp);
        [~, bestImagesEIG(w, r)] = max(EIG_temp);
        [~, bestImagesFT2(w, r)] = max(FT2_temp);

        % --------------------------------------------------------
        % Evaluate unimodality and monotonicity.
        % --------------------------------------------------------
        [unimodalVAR(w,r), monoVAR(w,r)] = ...
            evaluateFocusCurve(VAR_temp);

        [unimodalEIG(w,r), monoEIG(w,r)] = ...
            evaluateFocusCurve(EIG_temp);

        [unimodalFT2(w,r), monoFT2(w,r)] = ...
            evaluateFocusCurve(FT2_temp);

    end
end

%% Task 2 - Performance curves

figure;

plotIndex = 1;

for w = 1:numSizes

    for r = 1:numPositions

        subplot(numSizes, numPositions, plotIndex);

        VAR_curve = resultsVAR{w,r};
        EIG_curve = resultsEIG{w,r};
        FT2_curve = resultsFT2{w,r};

        plot(1:N, VAR_curve ./ max(VAR_curve), 'o-', 'LineWidth', 1.0);
        hold on;
        plot(1:N, EIG_curve ./ max(EIG_curve), 's-', 'LineWidth', 1.0);
        plot(1:N, FT2_curve ./ max(FT2_curve), '^-', 'LineWidth', 1.0);

        xlabel('Image');
        ylabel('Normalized measure');

        title(sprintf('%d x %d - %s', ...
              windowSizes(w), windowSizes(w), positionNames{r}));

        grid on;

        if plotIndex == 1
            legend('VAR', 'EIG', 'FT2', 'Location', 'best');
        end

        plotIndex = plotIndex + 1;

    end
end

sgtitle('Task 2 - Performance curves');

%% Task 2 - Best-focused images

fprintf('\n============================================================\n');
fprintf('TASK 2 - BEST FOCUS IMAGES\n');
fprintf('============================================================\n');

for w = 1:numSizes

    for r = 1:numPositions

        fprintf(['Window %d x %d, %s: VAR = %d, EIG = %d, ', ...
                 'FT2 = %d\n'], ...
                windowSizes(w), windowSizes(w), positionNames{r}, ...
                bestImagesVAR(w,r), ...
                bestImagesEIG(w,r), ...
                bestImagesFT2(w,r));

    end
end

%% Task 2 - Unimodality

fprintf('\n============================================================\n');
fprintf('TASK 2 - UNIMODALITY\n');
fprintf('============================================================\n');

fprintf(['A curve is classified as unimodal here if it is monotonically ', ...
         'increasing up to its global maximum and monotonically ', ...
         'decreasing after it.\n\n']);

for w = 1:numSizes

    for r = 1:numPositions

        fprintf('Window %d x %d, %s: ', ...
                windowSizes(w), windowSizes(w), positionNames{r});

        fprintf('VAR = %s, ', yesNo(unimodalVAR(w,r)));
        fprintf('EIG = %s, ', yesNo(unimodalEIG(w,r)));
        fprintf('FT2 = %s\n', yesNo(unimodalFT2(w,r)));

    end
end

%% Task 2 - Monotonicity

fprintf('\n============================================================\n');
fprintf('TASK 2 - MONOTONICITY\n');
fprintf('============================================================\n');

fprintf(['Monotonicity score: 1 means perfectly monotonic on both ', ...
         'sides of the peak; lower values indicate more deviations.\n\n']);

for w = 1:numSizes

    for r = 1:numPositions

        fprintf(['Window %d x %d, %s: VAR = %.3f, ', ...
                 'EIG = %.3f, FT2 = %.3f\n'], ...
                windowSizes(w), windowSizes(w), positionNames{r}, ...
                monoVAR(w,r), monoEIG(w,r), monoFT2(w,r));

    end
end

%% Task 2 - Computational efficiency

meanVARtime = mean(timeVAR(:));
meanEIGtime = mean(timeEIG(:));
meanFT2time = mean(timeFT2(:));

fprintf('\n============================================================\n');
fprintf('TASK 2 - COMPUTATIONAL EFFICIENCY\n');
fprintf('============================================================\n');

fprintf('Mean calculation time over all window tests:\n');
fprintf('VAR = %.6f s\n', meanVARtime);
fprintf('EIG = %.6f s\n', meanEIGtime);
fprintf('FT2 = %.6f s\n', meanFT2time);

% Display timing for each window size and position.
fprintf('\nCalculation times for each test:\n');

for w = 1:numSizes

    for r = 1:numPositions

        fprintf(['Window %d x %d, %s: VAR = %.6f s, ', ...
                 'EIG = %.6f s, FT2 = %.6f s\n'], ...
                windowSizes(w), windowSizes(w), positionNames{r}, ...
                timeVAR(w,r), timeEIG(w,r), timeFT2(w,r));

    end
end

figure;

bar([meanVARtime, meanEIGtime, meanFT2time]);

set(gca, 'XTickLabel', {'VAR', 'EIG', 'FT2'});
ylabel('Mean calculation time [s]');
title('Task 2 - Computational efficiency');
grid on;

%% Task 2 - Monotonicity comparison

meanMonoVAR = mean(monoVAR(:));
meanMonoEIG = mean(monoEIG(:));
meanMonoFT2 = mean(monoFT2(:));

figure;

bar([meanMonoVAR, meanMonoEIG, meanMonoFT2]);

set(gca, 'XTickLabel', {'VAR', 'EIG', 'FT2'});
ylabel('Mean monotonicity score');
ylim([0 1]);
title('Task 2 - Mean monotonicity');
grid on;

%% Task 2 - Summary and recommendation

meanUnimodalVAR = mean(unimodalVAR(:));
meanUnimodalEIG = mean(unimodalEIG(:));
meanUnimodalFT2 = mean(unimodalFT2(:));

meanBestVAR = mean(bestImagesVAR(:));
meanBestEIG = mean(bestImagesEIG(:));
meanBestFT2 = mean(bestImagesFT2(:));

fprintf('\n============================================================\n');
fprintf('TASK 2 - SUMMARY\n');
fprintf('============================================================\n');

fprintf('Mean unimodality rate:\n');
fprintf('VAR = %.3f\n', meanUnimodalVAR);
fprintf('EIG = %.3f\n', meanUnimodalEIG);
fprintf('FT2 = %.3f\n', meanUnimodalFT2);

fprintf('\nMean monotonicity:\n');
fprintf('VAR = %.3f\n', meanMonoVAR);
fprintf('EIG = %.3f\n', meanMonoEIG);
fprintf('FT2 = %.3f\n', meanMonoFT2);

fprintf('\nMean selected best image number:\n');
fprintf('VAR = %.2f\n', meanBestVAR);
fprintf('EIG = %.2f\n', meanBestEIG);
fprintf('FT2 = %.2f\n', meanBestFT2);

% A simple overall score combines the two requested curve criteria.
% This is used only to provide an objective summary of the experiments.
overallVAR = (meanUnimodalVAR + meanMonoVAR) / 2;
overallEIG = (meanUnimodalEIG + meanMonoEIG) / 2;
overallFT2 = (meanUnimodalFT2 + meanMonoFT2) / 2;

fprintf('\nOverall curve-quality score:\n');
fprintf('VAR = %.3f\n', overallVAR);
fprintf('EIG = %.3f\n', overallEIG);
fprintf('FT2 = %.3f\n', overallFT2);

scores = [overallVAR, overallEIG, overallFT2];
[~, recommendedMethod] = max(scores);

methodNames = {'VAR', 'EIG', 'FT2'};

fprintf('\nRecommendation based on unimodality + monotonicity:\n');
fprintf('%s\n', methodNames{recommendedMethod});

% Also report the fastest method.
times = [meanVARtime, meanEIGtime, meanFT2time];
[~, fastestMethod] = min(times);

fprintf('Fastest method in these experiments: %s\n', ...
        methodNames{fastestMethod});

fprintf('\nInterpretation:\n');
fprintf(['The performance curves should be inspected together with ', ...
         'the numerical scores above. A larger focus window generally ', ...
         'contains more image information but requires more computation. ', ...
         'The position also matters because different parts of the ', ...
         'scene can contain different amounts of texture/detail.\n']);

%% ================================================================
%  TASK 3 - CREATE A FULL FOCUS IMAGE
%  ================================================================

% Use VAR for the full-focus reconstruction, as in the original
% implementation.
windowSize = 32;
step = windowSize / 2;

tic;

[fullFocus, selectionMap] = createFullFocusImage( ...
    images, N, height, width, windowSize, step);

task3BaseTime = toc;

fprintf('\n============================================================\n');
fprintf('TASK 3 - FULL FOCUS IMAGE\n');
fprintf('============================================================\n');

fprintf('Window size: %d x %d\n', windowSize, windowSize);
fprintf('Step size: %d pixels\n', step);
fprintf('Calculation time: %.6f s\n', task3BaseTime);

figure;

imshow(fullFocus, []);
title(sprintf('Task 3 - Full-focus image using VAR (%d x %d)', ...
      windowSize, windowSize));

%% Task 3 - Compare different focus-window sizes

task3WindowSizes = [16 32 64];

task3Times = zeros(size(task3WindowSizes));
task3Results = cell(size(task3WindowSizes));
task3SelectionMaps = cell(size(task3WindowSizes));

for w = 1:length(task3WindowSizes)

    windowSize = task3WindowSizes(w);
    step = windowSize / 2;

    tic;

    [testFocus, testSelectionMap] = createFullFocusImage( ...
        images, N, height, width, windowSize, step);

    task3Times(w) = toc;

    task3Results{w} = testFocus;
    task3SelectionMaps{w} = testSelectionMap;

end

%% Task 3 - Display full-focus images

figure;

for w = 1:length(task3WindowSizes)

    subplot(1, length(task3WindowSizes), w);

    imshow(task3Results{w}, []);

    title(sprintf('%d x %d', ...
          task3WindowSizes(w), task3WindowSizes(w)));

end

sgtitle('Task 3 - Full-focus images for different window sizes');

%% Task 3 - Display selected image maps

% Each pixel shows which image in the focal stack was selected most
% often for the corresponding local focus window.
figure;

for w = 1:length(task3WindowSizes)

    subplot(1, length(task3WindowSizes), w);

    imagesc(task3SelectionMaps{w});
    axis image;
    colorbar;

    xlabel('x [pixels]');
    ylabel('y [pixels]');

    title(sprintf('Selected image map: %d x %d', ...
          task3WindowSizes(w), task3WindowSizes(w)));

end

sgtitle('Task 3 - Local best-focus image selection');

%% Task 3 - Computational time

fprintf('\n============================================================\n');
fprintf('TASK 3 - WINDOW SIZE AND COMPUTATION TIME\n');
fprintf('============================================================\n');

for w = 1:length(task3WindowSizes)

    fprintf('Window %d x %d: %.6f s\n', ...
            task3WindowSizes(w), ...
            task3WindowSizes(w), ...
            task3Times(w));

end

figure;

plot(task3WindowSizes, task3Times, 'o-', 'LineWidth', 1.5);

xlabel('Focus-window size [pixels]');
ylabel('Calculation time [s]');
title('Task 3 - Calculation time for different window sizes');
grid on;

%% Task 3 - Summary

fprintf('\nTask 3 interpretation:\n');
fprintf(['The smaller focus windows provide more local adaptation because ', ...
         'the best-focused image can change over shorter spatial distances. ', ...
         'However, smaller windows require more window positions and therefore ', ...
         'can increase the total calculation time. Larger windows require ', ...
         'fewer positions but can mix regions with different focus depths.\n']);

fprintf('\n============================================================\n');
fprintf('LAB 1 FINISHED\n');
fprintf('============================================================\n');


%% ================================================================
% LOCAL FUNCTIONS
% ================================================================

function [isUnimodal, monotonicity] = evaluateFocusCurve(curve)

    % Find the global maximum.
    [~, peak] = max(curve);

    % Values before and after the peak.
    left = curve(1:peak);
    right = curve(peak:end);

    % Handle the special case where the peak is at an endpoint.
    if length(left) > 1
        leftDiff = diff(left);
        leftScore = mean(leftDiff >= 0);
    else
        leftScore = 1;
    end

    if length(right) > 1
        rightDiff = diff(right);
        rightScore = mean(rightDiff <= 0);
    else
        rightScore = 1;
    end

    % Monotonicity score in [0,1].
    monotonicity = (leftScore + rightScore) / 2;

    % Strict unimodality criterion used for the binary result.
    if length(left) > 1
        leftMonotonic = all(diff(left) >= 0);
    else
        leftMonotonic = true;
    end

    if length(right) > 1
        rightMonotonic = all(diff(right) <= 0);
    else
        rightMonotonic = true;
    end

    isUnimodal = leftMonotonic && rightMonotonic;

end


function output = yesNo(value)

    if value
        output = 'Yes';
    else
        output = 'No';
    end

end


function [fullFocus, selectionMap] = createFullFocusImage( ...
    images, N, height, width, windowSize, step)

    fullFocus = zeros(height, width);
    weight = zeros(height, width);

    % Store the selected best image at every window position.
    selectionMap = zeros(height, width);
    selectionWeight = zeros(height, width);

    % Use a fixed window whenever possible.
    yPositions = 1:step:(height - windowSize + 1);
    xPositions = 1:step:(width - windowSize + 1);

    % Add the final possible window position if it is not already included.
    if yPositions(end) ~= height - windowSize + 1
        yPositions = [yPositions, height - windowSize + 1];
    end

    if xPositions(end) ~= width - windowSize + 1
        xPositions = [xPositions, width - windowSize + 1];
    end

    for y = yPositions

        for x = xPositions

            y2 = y + windowSize - 1;
            x2 = x + windowSize - 1;

            focusMeasure = zeros(1, N);

            % Evaluate VAR for every image in the focal stack.
            for i = 1:N

                window = images{i}(y:y2, x:x2);
                focusMeasure(i) = var(window(:));

            end

            % Select the image with the highest local variance.
            [~, bestImage] = max(focusMeasure);

            % Add the selected image region to the reconstruction.
            fullFocus(y:y2, x:x2) = ...
                fullFocus(y:y2, x:x2) + ...
                images{bestImage}(y:y2, x:x2);

            weight(y:y2, x:x2) = ...
                weight(y:y2, x:x2) + 1;

            % Store which focal-stack image was selected.
            selectionMap(y:y2, x:x2) = ...
                selectionMap(y:y2, x:x2) + bestImage;

            selectionWeight(y:y2, x:x2) = ...
                selectionWeight(y:y2, x:x2) + 1;

        end
    end

    % Average overlapping windows.
    fullFocus = fullFocus ./ weight;

    % Average the selected image numbers in overlapping windows.
    selectionMap = selectionMap ./ selectionWeight;

end

clear;
close all;
clc;

%% Load images
files = dir('KaktusFocus/*.jpg');

for i = 1:length(files)
    I = imread(fullfile('KaktusFocus',files(i).name));

    if size(I,3) == 3
        I = rgb2gray(I);
    end

    images{i} = double(I);
end

N = length(images);

%% Select focus window
figure;
imshow(images{13},[]);
title('Select two corners of the focus window');

[x,y] = ginput(2);

x1 = round(min(x));
x2 = round(max(x));
y1 = round(min(y));
y2 = round(max(y));

%% Calculate VAR, EIG and FT2
VAR = zeros(1,N);
EIG = zeros(1,N);
FT2 = zeros(1,N);

for i = 1:N

    I = images{i}(y1:y2,x1:x2);

    %% VAR
    VAR(i) = var(I(:));

    %% EIG
    Gx = imfilter(I,[-1 0 1],'replicate');
    Gy = imfilter(I,[-1;0;1],'replicate');

    EIG(i) = sum(Gx(:).^2 + Gy(:).^2);

    %% FT2
    F = fftshift(fft2(I));
    F = abs(F);

    [h,w] = size(F);
    [X,Y] = meshgrid(1:w,1:h);

    centerX = w/2;
    centerY = h/2;

    radius = sqrt((X-centerX).^2 + (Y-centerY).^2);

    FT2(i) = sum(F(radius > 10).^2);

end

%% Find best focus
[~,bestVAR] = max(VAR);
[~,bestEIG] = max(EIG);
[~,bestFT2] = max(FT2);

fprintf('Best VAR image = %d\n',bestVAR);
fprintf('Best EIG image = %d\n',bestEIG);
fprintf('Best FT2 image = %d\n',bestFT2);

%% Display best images
figure;

subplot(1,3,1);
imshow(images{bestVAR},[]);
title(['VAR: Image ',num2str(bestVAR)]);

subplot(1,3,2);
imshow(images{bestEIG},[]);
title(['EIG: Image ',num2str(bestEIG)]);

subplot(1,3,3);
imshow(images{bestFT2},[]);
title(['FT2: Image ',num2str(bestFT2)]);

%% Performance curves
figure;

plot(1:N,VAR/max(VAR),'o-');
hold on;
plot(1:N,EIG/max(EIG),'s-');
plot(1:N,FT2/max(FT2),'^-');

xlabel('Image number');
ylabel('Normalized focus measure');
legend('VAR','EIG','FT2');
title('Performance of focus measures');
grid on;


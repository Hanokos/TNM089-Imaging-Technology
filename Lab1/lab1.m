clear;
close all;
clc;
addpath(genpath('C:\TNM089\Project'))
%% Ladda bilder

files = dir('KaktusFocus/*.jpg');

for i = 1:length(files)
    I = imread(fullfile('KaktusFocus', files(i).name));

    if size(I,3) == 3
        I = rgb2gray(I);
    end

    images{i} = double(I);
end

N = length(images);

%% Välj fokusfönster

figure;
imshow(images{13},[]);
title('Välj två hörn för fokusfönstret');

[x,y] = ginput(2);

x1 = round(min(x));
x2 = round(max(x));
y1 = round(min(y));
y2 = round(max(y));

%% Beräkna VAR, EIG och FT2

VAR = zeros(1,N);
EIG = zeros(1,N);
FT2 = zeros(1,N);

for i = 1:N

    I = images{i}(y1:y2,x1:x2);

    % Variansen används som fokusmått
    VAR(i) = var(I(:));

    % Gradient används som fokusmått
    Gx = imfilter(I,[-1 0 1],'replicate');
    Gy = imfilter(I,[-1;0;1],'replicate');

    EIG(i) = sum(Gx(:).^2 + Gy(:).^2);

    % Fouriertransform används som fokusmått
    F = fft2(I);
    FT2(i) = sum(abs(F(:)));

end

%% Hitta bästa fokus

[~,bestVAR] = max(VAR);
[~,bestEIG] = max(EIG);
[~,bestFT2] = max(FT2);

fprintf('Bästa VAR-bild = %d\n',bestVAR);
fprintf('Bästa EIG-bild = %d\n',bestEIG);
fprintf('Bästa FT2-bild = %d\n',bestFT2);

%% Visa bästa bilder

figure;

subplot(1,3,1);
imshow(images{bestVAR},[]);
title(['VAR: Bild ',num2str(bestVAR)]);

subplot(1,3,2);
imshow(images{bestEIG},[]);
title(['EIG: Bild ',num2str(bestEIG)]);

subplot(1,3,3);
imshow(images{bestFT2},[]);
title(['FT2: Bild ',num2str(bestFT2)]);

%% Prestandakurvor

figure;

plot(1:N,VAR/max(VAR),'o-');
hold on;
plot(1:N,EIG/max(EIG),'s-');
plot(1:N,FT2/max(FT2),'^-');

xlabel('Bildnummer');
ylabel('Normaliserat fokusmått');
legend('VAR','EIG','FT2');
title('Prestanda för fokusmåtten');
grid on;

%% Task 3

windowSize = 32;
step = windowSize;

[height,width] = size(images{1});

fullFocus = zeros(height,width);

for y = 1:step:height-windowSize+1

    for x = 1:step:width-windowSize+1

        % Beräkna VAR för varje bild
        focusMeasure = zeros(1,N);

        for i = 1:N

            window = images{i}( ...
                y:y+windowSize-1, ...
                x:x+windowSize-1);

            focusMeasure(i) = var(window(:));

        end

        % Välj bilden med bäst fokus
        [~,bestImage] = max(focusMeasure);

        % Kopiera det bästa området
        fullFocus( ...
            y:y+windowSize-1, ...
            x:x+windowSize-1) = images{bestImage}( ...
            y:y+windowSize-1, ...
            x:x+windowSize-1);

    end
end

%% Visa fullfokusbilden

figure;
imshow(fullFocus,[]);
title('Fullfokusbild med VAR');


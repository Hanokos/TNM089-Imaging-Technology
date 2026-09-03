clear;
close all;
clc;

%% DEL 0 LADDA IN BILDERNA
% Läser in alla jpg-bilder i mappen KaktusFocus. Bilderna är redan
% gråskala, så vi gör bara om dem till double så vi kan räkna på dem.
files = dir('KaktusFocus/*.jpg');
images = cell(1,length(files));
for i = 1:length(files)
    I = imread(fullfile('KaktusFocus',files(i).name));
    images{i} = double(I);
end
N = length(images);
[H,W] = size(images{1});

%% DEL 1 VÄLJ FOKUSFÖNSTER
% Vi visar en av bilderna i stacken (bild 13 funkar bra som referens)
% och låter användaren klicka ut fönstret man vill titta på.
refIdx = 13;   % vilken bild vi använder för att välja fönster
% Rita ut en röda ram över valt fönster så man ser att det blev rätt
figure;
imshow(images{refIdx},[]);
title('Klicka i två hörn för att välja fokusfönster');
[x,y] = ginput(2);
x1 = round(min(x));
x2 = round(max(x));
y1 = round(min(y));
y2 = round(max(y));

% Se till att man inte klickar utanför bilden
x1 = max(1,x1); y1 = max(1,y1);
x2 = min(W,x2); y2 = min(H,y2);

%% Visa de bilder som respektive mått pekar ut som bäst
figure;
imshow(images{refIdx},[]);
hold on;
rectangle('Position',[x1,y1,x2-x1,y2-y1],'EdgeColor','r','LineWidth',2);
title('Valt fokusfönster');
hold off;

%% DEL 1/2 BERÄKNA FOKUSMÅTTEN (VAR, EIG, FT2) + TIDTAGNING PER MÅTT
VAR = zeros(1,N);
EIG = zeros(1,N);
FT2 = zeros(1,N);

% Vi tar tid på varje mått separat så vi kan jämföra hur snabba de är
tVAR = 0; tEIG = 0; tFT2 = 0;

for i = 1:N
    I = images{i}(y1:y2,x1:x2);   % klipper ut fokusfönstret ur bild i

    %% VAR bildvarians (Ekvation 1 i artikeln)
    % Enklaste måttet: hur mycket varierar pixelvärdena i fönstret?
    t0 = tic;
    VAR(i) = var(I(:));
    tVAR = tVAR + toc(t0);

    %% EIG energin i bildgradienten (Ekv. 4)
    % Tar fram kanter i x och yled och summerar kvadraten av dem.
    % Skarpa bilder har starkare kanter -> högre värde.
    t0 = tic;
    Gx = imfilter(I,[-1 0 1],'replicate');
    Gy = imfilter(I,[-1;0;1],'replicate');
    EIG(i) = sum(Gx(:).^2 + Gy(:).^2);
    tEIG = tEIG + toc(t0);

    %% FT2 spektrummått (Ekv. 8, enligt formeln)
    % Summerar hela Fourierspektrumets magnitud, förutom själva
    % DC termen (medelvärdet), precis som definitionen i artikeln säger.
    t0 = tic;
    F = abs(fft2(I));
    F(1,1) = 0;            % DC-termen ligger alltid i (1,1) innan fftshift
    FT2(i) = sum(F(:));    % summera resten av spektrumet
    tFT2 = tFT2 + toc(t0);
end

fprintf('--- Tid per fokusmått (%d bilder, fönster %dx%d) ---\n', ...
    N, y2-y1+1, x2-x1+1);
fprintf('VAR total tid: %.4f s (%.5f s/bild)\n', tVAR, tVAR/N);
fprintf('EIG total tid: %.4f s (%.5f s/bild)\n', tEIG, tEIG/N);
fprintf('FT2 total tid: %.4f s (%.5f s/bild)\n\n', tFT2, tFT2/N);

%% Hitta bästa fokus enligt varje mått
[~,bestVAR] = max(VAR);
[~,bestEIG] = max(EIG);
[~,bestFT2] = max(FT2);
fprintf('Bästa bild enligt VAR = %d\n',bestVAR);
fprintf('Bästa bild enligt EIG = %d\n',bestEIG);
fprintf('Bästa bild enligt FT2 = %d\n',bestFT2);

%% Prestandakurvor (normaliserade så alla toppar hamnar på 1)
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

figure;
plot(1:N,VAR/max(VAR),'o-');
hold on;
plot(1:N,EIG/max(EIG),'s-');
plot(1:N,FT2/max(FT2),'^-');
xlabel('Bildnummer');
ylabel('Normaliserat fokusmått');
legend('VAR','EIG','FT2');
title(sprintf('Fokusmåttens prestanda (fönster: [%d %d]-[%d %d])',x1,y1,x2,y2));
grid on;

%% DEL 2b VARIERA FÖNSTERSTORLEK OCH POSITION
%  (svarar på: "hur påverkar storlek och position resultatet?")
% Utgångspunkt: samma fönster som valdes ovan med ginput
baseW = x2 - x1 + 1;
baseH = y2 - y1 + 1;
cx = round((x1+x2)/2);
cy = round((y1+y2)/2);

% Effekten av FÖNSTERSTORLEK (samma mittpunkt hela tiden)
sizeFactors = [0.5 1 2 4];   % hur mycket vi skalar upp/ner ursprungsfönstret
figure('Name','Effekt av fönsterstorlek');
for s = 1:length(sizeFactors)
    wS = round(baseW*sizeFactors(s));
    hS = round(baseH*sizeFactors(s));
    xa = max(1, cx-round(wS/2)); xb = min(W, xa+wS-1);
    ya = max(1, cy-round(hS/2)); yb = min(H, ya+hS-1);

    v = zeros(1,N); e = zeros(1,N); f = zeros(1,N);
    for i = 1:N
        I = images{i}(ya:yb,xa:xb);
        v(i) = var(I(:));
        Gx = imfilter(I,[-1 0 1],'replicate');
        Gy = imfilter(I,[-1;0;1],'replicate');
        e(i) = sum(Gx(:).^2 + Gy(:).^2);
        Fm = abs(fft2(I));
        Fm(1,1) = 0;
        f(i) = sum(Fm(:));
    end

    subplot(2,2,s);
    plot(1:N, v/max(v), 'o-'); hold on;
    plot(1:N, e/max(e), 's-');
    plot(1:N, f/max(f), '^-');
    title(sprintf('%dx%d fönster (faktor %.1fx)', yb-ya+1, xb-xa+1, sizeFactors(s)));
    xlabel('Bildnummer'); ylabel('Normaliserat mått');
    legend('VAR','EIG','FT2','Location','southoutside','Orientation','horizontal');
    grid on;
end
sgtitle('Hur fönsterstorleken påverkar fokusmåttens kurvor');

% Effekten av FÖNSTERPOSITION (samma storlek, olika platser i bilden)
% OBS: "center" motsvarar vårt egna valda fönster. De andra tre är
% helt nya positioner (fasta procentandelar av bildens bredd/höjd),
% inte förskjutningar av ditt fönster.
positions = struct( ...
    'name', {'centrum','uppe-vänster','uppe-höger','nere-mitten'}, ...
    'cx',   {cx, round(W*0.2), round(W*0.8), round(W*0.5)}, ...
    'cy',   {cy, round(H*0.2), round(H*0.2), round(H*0.8)});

figure('Name','Effekt av fönsterposition');
for p = 1:length(positions)
    xa = max(1, positions(p).cx-round(baseW/2)); xb = min(W, xa+baseW-1);
    ya = max(1, positions(p).cy-round(baseH/2)); yb = min(H, ya+baseH-1);

    v = zeros(1,N); e = zeros(1,N); f = zeros(1,N);
    for i = 1:N
        I = images{i}(ya:yb,xa:xb);
        v(i) = var(I(:));
        Gx = imfilter(I,[-1 0 1],'replicate');
        Gy = imfilter(I,[-1;0;1],'replicate');
        e(i) = sum(Gx(:).^2 + Gy(:).^2);
        Fm = abs(fft2(I));
        Fm(1,1) = 0;
        f(i) = sum(Fm(:));
    end

    subplot(2,2,p);
    plot(1:N, v/max(v), 'o-'); hold on;
    plot(1:N, e/max(e), 's-');
    plot(1:N, f/max(f), '^-');
    title(positions(p).name);
    xlabel('Bildnummer'); ylabel('Normaliserat mått');
    legend('VAR','EIG','FT2','Location','southoutside','Orientation','horizontal');
    grid on;
end
sgtitle('Hur fönsterpositionen påverkar fokusmåttens kurvor');


%%  DEL 3 FULLFOKUSBILD (sätter ihop bilden via ett litet glidande fönster)
%  Uppgift 3: "loopa över bilden med ett litet fokusfönster och använd,
%  för varje position, bilddatan från den stackbild som har bäst fokus."
%  Vi visar hur fönsterstorleken påverkar både bildkvalitet
%  och beräkningstid.
%%  3a. En fullfokusbild med en standardstorlek på fönstret
method       = 'VAR';   % 'VAR' | 'EIG' | 'FT2'
windowHeight = 16;
windowWidth  = 16;

[focusedImage, time] = buildFullFocusImage(images, N, H, W, ...
    windowHeight, windowWidth, method);

figure;
imshow(focusedImage,[]);
title(['Fullfokusbild - ' method ' (' num2str(windowHeight) 'x' num2str(windowWidth) ')']);
fprintf('\nMetod: %s\n',method);
fprintf('Fönster: %d x %d\n',windowHeight,windowWidth);
fprintf('Tid: %.3f sekunder\n',time);

%% 3b. Variera fönsterstorlek: kvalitet vs beräkningstid
windowSizes = [4 8 16 32 64];   % kvadratiska fönster, från litet till stort
sweepTimes  = zeros(1,length(windowSizes));
sweepImages = cell(1,length(windowSizes));

for wsI = 1:length(windowSizes)
    ws = windowSizes(wsI);
    [sweepImages{wsI}, sweepTimes(wsI)] = buildFullFocusImage(images, N, H, W, ...
        ws, ws, method);
    fprintf('Fönster %3dx%-3d -> tid = %.3f s\n', ws, ws, sweepTimes(wsI));
end

% Visa alla fönsterstorlekar sida vid sida så man kan jämföra kvaliteten
figure('Name','Fullfokusbild vid olika fönsterstorlekar');
nCols = length(windowSizes);
for wsI = 1:nCols
    subplot(1,nCols,wsI);
    imshow(sweepImages{wsI},[]);
    title(sprintf('%dx%d\n%.2fs', windowSizes(wsI), windowSizes(wsI), sweepTimes(wsI)));
end
sgtitle(sprintf('Kvalitet vs fönsterstorlek (metod: %s)', method));

% Plotta hur beräkningstiden beror på fönsterstorleken
figure('Name','Beräkningstid vs fönsterstorlek');
plot(windowSizes, sweepTimes, 'o-','LineWidth',1.5);
xlabel('Fönsterstorlek (pixlar, kvadratiskt)');
ylabel('Beräkningstid (s)');
title(sprintf('Beräkningstid vs fönsterstorlek, del 3 (metod: %s)', method));
grid on;

%% LOKAL FUNKTION - bygger en fullfokusbild för en given fönsterstorlek
function [focusedImage, elapsed] = buildFullFocusImage(images, N, H, W, ...
    windowHeight, windowWidth, method)

    focusedImage = zeros(H,W);

    % Startpositioner för raderna/kolumnerna. Vi ser till att sista
    % blocket alltid når hela vägen till kanten, även om H eller W
    % inte är jämnt delbart med fönsterstorleken.
    yStarts = 1:windowHeight:H;
    xStarts = 1:windowWidth:W;

    tic
    for yi = 1:length(yStarts)
        y = yStarts(yi);
        yEnd = min(y+windowHeight-1, H);   % klipp vid nedre kanten

        for xi = 1:length(xStarts)
            x = xStarts(xi);
            xEnd = min(x+windowWidth-1, W); % klipp vid högra kanten

            % Testa alla N bilder på just den här blockpositionen
            % och se vilken som är skarpast där
            scores = zeros(1,N);
            for k = 1:N
                I = images{k}(y:yEnd, x:xEnd);

                if strcmp(method,'VAR')
                    scores(k) = var(I(:));
                elseif strcmp(method,'EIG')
                    Gx = imfilter(I,[-1 0 1],'replicate');
                    Gy = imfilter(I,[-1;0;1],'replicate');
                    scores(k) = sum(Gx(:).^2 + Gy(:).^2);
                elseif strcmp(method,'FT2')
                    Fm = abs(fft2(I));
                    Fm(1,1) = 0;          % nollställ bara DC-termen (Ekv. 8)
                    scores(k) = sum(Fm(:));
                end
            end

            % Plocka blocket från den bild som fick högst poäng
            [~,bestImage] = max(scores);
            focusedImage(y:yEnd, x:xEnd) = images{bestImage}(y:yEnd, x:xEnd);
        end
    end
    elapsed = toc;
end

%% flashNoFlashRun.m -- compare NO-FLASH INPUT vs our BEST RESULT (denoised + detail transfer)
% Full image resolution, no resizing. Needs Image Processing Toolbox.
close all; clear; clc;
imgDir  = uigetdir(pwd, 'Select the folder with the flash/no-flash images');
if isequal(imgDir, 0), error('No folder selected.'); end

names       = {'lamp','cave01','carpet','puppets','potsdetail'};
stepThrough = true;    % press a key between scenes (figures close)
numRegions  = 2;       % zoom regions to pick per scene

%% ---- TUNING KNOBS (full-resolution pixels) ----
sdDenoise  = 6;      % spatial sigma for denoising
srJoint    = 0.02;   % joint bilateral range sigma (larger = smoother)
sdDetail   = 10;     % spatial sigma of detail layer (larger = more detail transferred)
srDetail   = 0.10;   % range sigma of detail layer   (larger = more detail, more halos)
detailGain = 1.3;    % 1 = paper; >1 exaggerates detail
epsv       = 0.02;
diffMax    = 0.15;   % change (0-1 scale) that shows as full white in the difference image

for n = 1:numel(names)
    nm = names{n};
    A = im2double(imread(fullfile(imgDir,[nm '_01_noflash.tif'])));   % no-flash input
    F = im2double(imread(fullfile(imgDir,[nm '_00_flash.tif'])));     % flash (used as guide only)
    if size(A,3)==1, A = repmat(A,[1 1 3]); F = repmat(F,[1 1 3]); end
    A = min(max(A,0),1);  F = min(max(F,0),1);
    [H,W,~] = size(A);
    fprintf('\n%s: processing full size %d x %d ...\n', nm, W, H);

    sdD = sdDenoise; sjr = srJoint; sdT = sdDetail; srT = srDetail;
    if strcmp(nm,'cave01'), sdD = 12; sjr = 0.05; sdT = 16; srT = 0.15; end  % very noisy pair

    tic;
    Anr   = fastJbf(A, F, sdD, sjr);                               % joint bilateral denoise
    Fbase = fastJbf(F, F, sdT, srT);                               % smooth flash -> detail layer
    Afin  = min(max(Anr .* (((F+epsv)./(Fbase+epsv)).^detailGain), 0), 1);   % best result
    dif   = mean(abs(Afin - A), 3);                                % what changed
    fprintf('  done in %.1f s | mean abs change vs no-flash: %.4f\n', toc, mean(dif(:)));

    titles = {'INPUT: No-flash (noisy)', 'RESULT: Denoised + detail transfer', ...
              sprintf('DIFFERENCE (brighter = bigger change, mean %.3f)', mean(dif(:)))};

    % ---- overview (screen scaling only; data stays full resolution) ----
    figure('Name',[nm ' overview'],'NumberTitle','off','Color','w','Units','normalized','Position',[0.02 0.2 0.96 0.55]);
    for k = 1:3
        subplot(1,3,k);
        if k < 3, imshow(A*(k==1) + Afin*(k==2)); else, imshow(dif,[0 diffMax]); colormap(gca,'hot'); end
        title(titles{k},'FontSize',11);
    end
    sgtitle([nm ' - overview'],'FontWeight','bold');
    drawnow;

    % ---- zoom regions (full-resolution pixels), panels linked ----
    for z = 1:numRegions
        fsel = figure('Name','Pick region','NumberTitle','off');
        imshow(Afin);
        title(sprintf('%s: region %d/%d - DRAG a rectangle, then double-click inside it', nm, z, numRegions));
        rect = getrect(gca);  close(fsel);
        if isempty(rect) || rect(3) < 3 || rect(4) < 3
            rect = [W*0.4 H*0.4 W*0.15 H*0.15];
        end
        xl = [max(1,rect(1)), min(W,rect(1)+rect(3))];
        yl = [max(1,rect(2)), min(H,rect(2)+rect(4))];

        figure('Name',sprintf('%s zoom %d',nm,z),'NumberTitle','off','Color','w', ...
               'Units','normalized','Position',[0.02 0.2 0.96 0.55]);
        ax = gobjects(1,3);
        for k = 1:3
            ax(k) = subplot(1,3,k);
            if k == 1, imshow(A); elseif k == 2, imshow(Afin);
            else, imshow(dif,[0 diffMax]); colormap(ax(k),'hot'); end
            title(titles{k},'FontSize',11);
        end
        linkaxes(ax,'xy');
        xlim(ax(1),xl); ylim(ax(1),yl);
        sgtitle(sprintf('%s - region %d, full resolution (zoom/pan: all panels stay in sync)', nm, z),'FontWeight','bold');
        drawnow;
    end
    if stepThrough && n < numel(names)
        fprintf('Inspect the figures, then press any key in the Command Window for the next scene...\n');
        pause; close all;
    end
end

%% ---------------------------------------------------------------
function out = fastJbf(A, G, sd, sr)
% Fast joint bilateral filter, piecewise-linear approximation (Durand & Dorsey 2002).
% With G = A this is the basic bilateral filter.
A = single(A);  G = single(G);
nL  = min(101, max(6, ceil(1/sr) + 1));      % level spacing <= sr
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
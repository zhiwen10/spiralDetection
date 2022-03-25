githubdir2 = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir2, 'spikes')))
addpath(genpath(fullfile(githubdir2, 'Pipelines')))
addpath(genpath(fullfile(githubdir2, 'NeuroPattToolbox'))) %https://github.com/BrainDynamicsUSYD/NeuroPattToolbox
addpath(genpath(fullfile(githubdir2, 'widefield')))
addpath(genpath(fullfile(githubdir2, 'npy-matlab')))
%% SVD, plot example trace, overlay with pupil and stim time
% specify folder for widefield image data
mn = 'ZYE_0012';
td = '2020-10-16';
en = 5;
serverRoot = expPath(mn, td, en);
wf_svd;
%%
downscale = 4; % downscale factor for each frame to process faster, by sacrificing resolution
lowpass = 0; % 1 is low pass filter <0.2Hz, 0 is bandpass filter between 2-8Hz
tStart = 1000; tEnd = 2000; % find spirals between time tStart:tEnd
useGPU = 1; % useGPU to process HS flowfield
Fs = 35; % frame sampling rate
params = setNeuroPattParams_mod(Fs); % params for NeuroPattToolbox
params.downsampleScale = 1; % no downsample in NeuroPattToolbox 
params.minCritRadius = 8; % spirals has to be at least this value of pixels in radius
checkRadius = 2:8; % spirals has to be at least this value of pixels in radius
%%
U1 = U(1:downscale:end,1:downscale:end,1:50);
dV1 = dV(1:50,tStart:tEnd);
t1 = t(tStart:tEnd);
%%
[trace2d,traceAmp,tracePhase] = spiralPhaseMap(U1,dV1,t1,lowpass);
%% check example trace
% figure; 
% plot(t,trace2d(:,40,40))
%% apply mask, this helps speed up spiral detection later
mimg1 = mimg(1:downscale:end,1:downscale:end);
figure; imagesc(mimg1);
BW = roipoly;
%%
mimg1(not(BW)) = 0;
figure; imagesc(mimg1);
%% HS flowfield calculation with GPU
[vxRaw,vyRaw] = HS_flowfield(tracePhase,useGPU);
%% example neighbour frames and flow field map for sanity check 
figure;
frame = 500;
AOtRaw1 = squeeze(tracePhase(frame,:,:));
AOtRaw2 = squeeze(tracePhase(frame+1,:,:));
vxRaw2 = zeros(size(vxRaw,2,3));
vyRaw2 = zeros(size(vyRaw,2,3));
skip = 2;
vxRaw2(1:skip:end,1:skip:end) = vxRaw(frame,1:skip:end,1:skip:end);
vyRaw2(1:skip:end,1:skip:end) = vyRaw(frame,1:skip:end,1:skip:end);
vscale = 2;
subplot(1,3,1)
imHORaw1 = imagesc(AOtRaw1);
axis image; axis off;
colormap(hsv);
subplot(1,3,2)
imHORaw2 = imagesc(AOtRaw2);
axis image; axis off;
colormap(hsv);
subplot(1,3,3)
imHORaw3 = imagesc(AOtRaw1);
axis image; axis off;
colormap(hsv);
hold on;
imH1Raw4 = quiver(vxRaw2*vscale,vyRaw2*vscale,'k','lineWidth',1,'autoScale','off');
%% spiral detection preparation
vxRaw = permute(vxRaw,[2,3,1]);
vyRaw = permute(vyRaw,[2,3,1]);
thisvfRaw = vxRaw + 1i*vyRaw; 
thisvfRaw = reshape(thisvfRaw,size(thisvfRaw,1)*size(thisvfRaw,2), size(thisvfRaw,3));
% mask out the area outside of the brain speed up the process of finding cirtical points
thisvfRaw(not(BW(:)),:) = nan;
x = size(vxRaw,1); y = size(vxRaw,2);
%% spiral detection
step = 500;
Nstep = ceil((length(t1)-1)/step); % flowfield vectors are one frame less
pattLocs = [];
tic
for ibatch = 1:Nstep
    if ibatch <Nstep
        thisvfRaw1 = thisvfRaw(:,1+(ibatch-1)*step:step*ibatch);
    elseif ibatch== Nstep
        thisvfRaw1 = thisvfRaw(:,1+(ibatch-1)*step:end);
    end
    thisvfRaw1 = reshape(thisvfRaw1,x,y,size(thisvfRaw1,2));
    [pattTypesRaw,pattLocsRaw] = ...
        findAllPatterns2(real(thisvfRaw1), imag(thisvfRaw1), params, ...
        angle(thisvfRaw1));
    pattLocs = cat(1,pattLocs,pattLocsRaw);
    fprintf('Step %g/%g; time elapsed %g seconds \n', [ibatch, Nstep, toc])
end
spirals = cat(1,pattLocs{:});
%% check if within all checkRadius has the same phase change signs
LocsSpiral = spiralFirstCheck(pattLocs,tracePhase,checkRadius);
spirals1 = cat(1,LocsSpiral{:});
%% check if all [flow vector sum of half of the spiral] within checkRadius have all same signs
LocsFlow = spiralDoubleCheck(vxRaw,vyRaw,pattLocs,checkRadius);
spirals2 = cat(1,LocsFlow{:});
%% meet both above conditions, and have same spiral direction
[spiralsAll,ia,ib] = intersect(spirals1(:,[1:4,7,8]),spirals2(:,[1:4,7,8]),'rows');
spiralsAll1 = spirals1(ia,:);
%% 
color1 = zeros(size(spiralsAll1,1),3);
color1(:,2) = spiralsAll1(:,7);
color1(:,1) = not(spiralsAll1(:,7));
figure;
imagesc(mimg1)
axis off; axis image;
colormap(gray)
hold on;
scatter(spiralsAll1(:,2),spiralsAll1(:,1),3,color1,'filled');
set(gca, 'ydir','reverse')
%% reorganize spirals to cell array, ordered by frames
frameMax = max(spiralsAll1(:,8));
for i = 1:frameMax
    ispirals = find(spiralsAll1(:,8)==i);
    if not(isempty(ispirals))
        spiralsCell{i} = spiralsAll1(ispirals,:);
    else
        spiralsCell{i} = [];
    end
end
%% spiralVideo Example
frameN = 200;
tracePhase1 = tracePhase(1:frameN ,:,:);
vxRaw1 = zeros(size(vxRaw,1),size(vxRaw,2),frameN); vyRaw1 = zeros(size(vyRaw,1),size(vyRaw,2),frameN);
vxRaw1(1:4:end,1:4:end,:) = vxRaw(1:4:end,1:4:end,1:frameN); vyRaw1(1:4:end,1:4:end,:) = vyRaw(1:4:end,1:4:end,1:frameN);
spiralsCell1 = spiralsCell(1:frameN );
fname = 'testVideo';
frameTAll = t1(1:200);
vscale = 4;
spiralVideo(tracePhase1,spiralsCell1,vxRaw1,vyRaw1,vscale,frameTAll,fname);
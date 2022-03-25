%%
figure; 
% imagesc(mimg(1:8:end,1:8:end))
% colormap(gray)
% hold on;
scatter_kde(spirals2(:,2),spirals2(:,1),'filled','MarkerSize',4);
set(gca, 'ydir','reverse')
%%
roi = drawpolygon;
%%
tf = inROI(roi,spirals2(:,2),spirals2(:,1));
%%
spiralAmpMap = traceAmp(tf,:,:);
meanAmpMap = squeeze(mean(spiralAmpMap,1));
%%
figure; imagesc(meanAmpMap);
axis off; axis image;
%%
p = randperm(size(spiralAmpMap,1),100);
maxi = max(spiralAmpMap(:)); mini = min(spiralAmpMap(:));
figure; 
for i = 1:100
    subplottight(10,10,i)
    imagesc(squeeze(spiralAmpMap(p(i),:,:)));
    % caxis([mini,maxi]);
    axis off; axis image;
end
%%
spiralAmpMap1 = reshape(spiralAmpMap,size(spiralAmpMap,1),size(spiralAmpMap,2)*size(spiralAmpMap,3));
spiralAmpMap1= spiralAmpMap1';
[Us,S,Vs]= svd(spiralAmpMap1,'econ');
%%
S1 = diag(S);
figure;
scatter(1:length(S1),S1);
%%
maxi = max(Us(:)); mini = min(Us(:));
Us = reshape(Us,x,y,size(Us,2));
figure; 
for i =1:10
    subplottight(1,10,i)
    imagesc(Us(:,:,i))
    caxis([mini,maxi]);
    axis off; axis image;
end
%%
traceAmp = reshape(traceAmp,size(traceAmp,1),size(traceAmp,2)*size(traceAmp,3));
traceAmp= traceAmp';
[Us,S,Vs]= svd(traceAmp,'econ');
%%
S1 = diag(S);
figure;
scatter(1:length(S1),S1);
%%
Us = reshape(Us,x,y,size(Us,2));
%%
Us1 = Us(:,:,1:10);
maxi = max(Us1(:)); mini = min(Us1(:));
figure; 
for i =1:10
    subplottight(1,10,i)
    imagesc(Us1(:,:,i))
    caxis([mini,maxi]);
    axis off; axis image;
end
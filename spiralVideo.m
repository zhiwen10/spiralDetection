
function spiralVideo(tracePhase1,pattLocs,vx2Raw,vy2Raw,vscale,frameTAll,fname)
% vscale = 5;
frameN = length(pattLocs);
minO = min(tracePhase1(:));
maxO = max(tracePhase1(:));
ha = figure;
ha.Position = [100 100 360 360];
% fname = 'ZYE12-HS_phase_globalWithSpiral_dV3.avi';
v = VideoWriter(fname);
v.FrameRate = 7;
xlimV1 = floor(size(tracePhase1,2)*0.1);
xlimV = [-xlimV1;size(tracePhase1,2)+xlimV1];
ylimV1 = floor(size(tracePhase1,3)*0.1);
ylimV = [-xlimV1;size(tracePhase1,3)+ylimV1];
open(v);
for i = 1:frameN
    AOtRaw = squeeze(tracePhase1(i,:,:));
    ispiralRaw = pattLocs{i};
    if isempty(ispiralRaw)
        ispiralRaw = zeros(1,8);
    end
    frameT = frameTAll(i);
    if i == 1
        imHORaw = imagesc(AOtRaw);
        axis image; axis off;
        caxis([minO, maxO])
        colormap(hsv);
        hold on;
        imH1Raw = quiver(vx2Raw(:,:,i)*vscale,vy2Raw(:,:,i)*vscale,'k','lineWidth',1,'autoScale','off');
        hold on;
        imH2Raw= scatter(ispiralRaw(:,2),ispiralRaw(:,1),20,'k','filled');
        if ispiralRaw(1,1)==0
            set(imH2Raw,'MarkerFaceAlpha',0); set(imH2Raw,'MarkerEdgeAlpha',0);
        else
            set(imH2Raw,'MarkerFaceAlpha',1); set(imH2Raw,'MarkerEdgeAlpha',1);
        end
        set(gca, 'ydir','reverse')
        xlim(xlimV); ylim(ylimV);
        h1 = text(1,46,{['frame: ' num2str(i)],['Time: ' num2str(frameT)]},'FontSize',6,'color','white');
    else
        set(imHORaw, 'CData', AOtRaw);
        set(imH1Raw,'UData', squeeze(vx2Raw(:,:,i)*vscale));
        set(imH1Raw,'VData', squeeze(vy2Raw(:,:,i)*vscale));
        set(imH2Raw,'XData', ispiralRaw(:,2));
        set(imH2Raw,'YData', ispiralRaw(:,1));
        h1.String = {['frame: ' num2str(i)],['Time: ' num2str(frameT)]};
        if ispiralRaw(1,1)==0
            set(imH2Raw,'MarkerFaceAlpha',0); set(imH2Raw,'MarkerEdgeAlpha',0);
        else
            set(imH2Raw,'MarkerFaceAlpha',1); set(imH2Raw,'MarkerEdgeAlpha',1);
        end
    end   
    thisFrame = getframe(ha);
    writeVideo(v, thisFrame);
end
close(v)
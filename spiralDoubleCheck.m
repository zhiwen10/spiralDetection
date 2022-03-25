function LocsFlow = spiralDoubleCheck(vxRaw,vyRaw,pattLocs,checkRadius)
x = size(vxRaw,1); y = size(vxRaw,2);
frameN = length(pattLocs);
%%
tic
for frame = 1:frameN
    % column 3 and 4 of pattLocs are in and out spirals
    mpoint = cat(1,pattLocs{frame,3:4}); 
    vfx = squeeze(vxRaw(:,:,frame));   
    vfy = squeeze(vyRaw(:,:,frame));   
    if ~isempty(mpoint)
        for kk = 1:size(mpoint,1)
            num1 = 1;
            pwSign = zeros(1,6);
            pwRadius = zeros(1,6);
            for rs1 = checkRadius
                px3 = mpoint(kk,2); py3 = mpoint(kk,1); r3 = rs1;
                [pwsign] = spiralCheckFlow(vfx,vfy,px3,py3,x,y,r3);
                if isempty(pwsign)
                    r3 = nan;
                    pwsign = nan;
                end
                pwRadius(num1) = r3;
                pwSign(num1) = pwsign;
                num1 = num1+1;
            end
            % if sum of all flow direction of a half circle are positive, then it's conterclockwise
            if not(all(isnan(pwSign))) && all(pwSign>0)
                mpoint(kk,7) = 1;
            % if sum of all flow direction of a half circle are negative, then it's clockwise
            elseif not(all(isnan(pwSign))) && all(pwSign<0)
                mpoint(kk,7) = 0;
            % if not consistent, then it's not a spiral
            else
                mpoint(kk,7) = nan;
            end               
        end
        mpoint(isnan(mpoint(:,7)),:) = [];
        mpoint(:,8) = frame;
        LocsFlow{frame} = mpoint;
        
    else
        LocsFlow{frame} = [];
    end
    if not(mod(frame,1000))
        fprintf('frame %g/%g; time elapsed %g seconds \n', [frame, frameN, toc])
    end
end

%%
spirals1 = cat(1,LocsFlow{:});
%%
color1 = zeros(size(spirals1,1),3);
color1(:,2) = spirals1(:,7);
color1(:,1) = not(spirals1(:,7));
figure;
% imagesc(mimg(1:8:end,1:8:end))
% colormap(gray)
% hold on;
scatter(spirals1(:,2),spirals1(:,1),3,color1,'filled');
set(gca, 'ydir','reverse')
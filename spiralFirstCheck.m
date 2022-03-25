function LocsSpiral = spiralFirstCheck(pattLocs,tracePhase,checkRadius)
% all phase angle differences around a spiral circle should be same
% direction. 
%%
frameN = size(pattLocs,1);
tic
for frame = 1:size(pattLocs,1)
        % column 3 and 4 of pattLocs are in and out spirals
        mpoint = cat(1,pattLocs{frame,3:4});       
        A = squeeze(tracePhase(frame,:,:));    
        if ~isempty(mpoint)
            for kk = 1:size(mpoint,1)
                num1 = 1;
                for rs1 = checkRadius
                    px3 = mpoint(kk,2); py3 = mpoint(kk,1); r3 = rs1;
                    [pw,pwsign] = spiralCheckpw1(A,px3,py3,r3);
                    if isempty(pwsign)
                        r3 = nan;
                        pwsign = nan;
                    end
                    pwRadius(num1) = r3;
                    pwSign(num1) = pwsign;
                    num1 = num1+1;
                end
                mpoint(kk,5) = min(pwRadius); mpoint(kk,6) = max(pwRadius); 
                % if all signs are real, then it's conterclockwise spiral
                if not(all(isnan(pwSign))) && all(pwSign)
                    mpoint(kk,7) = 1;
                % if all signs are not real, then it's clockwise spiral
                elseif not(all(isnan(pwSign))) && not(any(pwSign))
                    mpoint(kk,7) = 0;
                % otherwise not a spiral
                else
                    mpoint(kk,7) = nan;
                end               
            end
            mpoint(isnan(mpoint(:,7)),:) = [];
            mpoint(:,8) = frame;
            LocsSpiral{frame} = mpoint;
        else
            LocsSpiral{frame} = [];
        end
        frame = frame+1;
        if not(mod(frame,100))
            fprintf('frame %g/%g; time elapsed %g seconds \n', [frame, frameN, toc])
        end
end

%%
locsSpiral1 = LocsSpiral;
spirals1 = cat(1,locsSpiral1{:});
%%
% size filter; only acccept spiral with x radius
spirals2 = spirals1;
% spirals2 = spirals1(spirals1(:,6)-spirals1(:,5)>=2,:);
color1 = zeros(size(spirals2 ,1),3);
color1(:,1) = spirals2(:,7);
color1(:,2) = double(not(spirals2(:,7)));
%%
figure;
% imagesc(mimg(1:8:end,1:8:end))
% colormap(gray)
% hold on; 
scatter(spirals2(:,2),spirals2(:,1),3,color1,'filled');
set(gca, 'ydir','reverse')
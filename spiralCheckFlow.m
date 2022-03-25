function [pwsignx] = spiralcheckFlow(vfx,vfy,px,py,x,y,r)
th = 30:30:150;
cy = round(r*sind(th)+py); 
cx = round(r*cosd(th)+px);
if all(cy>0) && all(cy<y) && all(cx>0) && all(cx<x)
    indx = sub2ind(size(vfx),cy, cx);
    phx = vfx(indx);
    phy = vfy(indx);
    pwsignx = sign(sum(phx));
    pwsigny = sign(sum(phy));
else
    pwsignx = nan;
end
if pwsignx>0
    color1 = 'r';
elseif pwsignx<0
    color1 = 'g';
else 
    color1 = 'k';
end
% scatter(cx,cy,10,color1,'filled');
% scatter(mpoint(kk,2),mpoint(kk,1),5,'b');
% quiver(cx,cy,phx,phy)
end
function [pw,pwsign] = spiralCheckpw1(A,px,py,r)
xsize = size(A,1); ysize = size(A,2);
th = 1:6:360; 
gw = myGaussWin(15, 1/6); 
cx = round(r*sind(th)+px); 
cy = round(r*cosd(th)+py);
if all(cy>0) && all(cy<ysize) && all(cx>0) && all(cx<xsize)
    ph = A(sub2ind(size(A),cy, cx));
    phdiff = angdiff(ph); 
    phdiff = conv(phdiff, gw, 'same'); 
    if all(phdiff>0)
        % red if a positive (counterclockwise)
        pw = [px, py]; pwsign = 0; 
    elseif all(phdiff<0)
        pw = [px, py]; pwsign = 1;     
    else
        pw = []; pwsign = [];    
    end
else
    pw = []; pwsign = []; 
end
end
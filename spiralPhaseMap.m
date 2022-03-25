function [trace2d,traceAmp,tracePhase] = spiralPhaseMap(U1,dV1,t,lowpass)
%%
nSV = size(U1,3);
Fs = 1/median(diff(t));
x = size(U1,1); y = size(U1,2);
%%
Ur = reshape(U1, size(U1,1)*size(U1,2), size(U1,3));
%% trace re-construction
meanTrace = Ur*dV1;
meanTrace = double(meanTrace);
tsize = size(meanTrace,2);
%% filter 2-8Hz
meanTrace = meanTrace -mean(meanTrace ,2);
% filter and hilbert transform work on each column
meanTrace = meanTrace';
if lowpass
    [f1,f2] = butter(6, 0.2/(Fs/2), 'low');
else
    [f1,f2] = butter(6, [2 8]/(Fs/2), 'bandpass');
end
meanTrace = filtfilt(f1,f2,meanTrace);
traceHilbert =hilbert(meanTrace);
tracePhase = angle(traceHilbert);
traceAmp = abs(traceHilbert);
%%
tracePhase = reshape(tracePhase,tsize,x,y);
traceAmp = reshape(traceAmp,tsize,x,y);
trace2d = reshape(meanTrace,tsize,x,y);
end
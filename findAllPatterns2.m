function [pattTypes,critpointLocs] = ...
    findAllPatterns2(vfx, vfy, params, phase, repNumber)
% FINDALLPATTERNS finds the patterns at each point in time for
% the velocity fields determined by the XxYxTime matrices VFX and
% VFY.
% PATTERNS is a matrix containing the type, location and time of each
%   pattern present: each pattern is a row, column titles are given by
%   COLNAMES. Pattern types are given by indices, which correspond to the
%   patterns named in PATTTYPES.
% ALLPATTERNLOCS is a Nx1 cell array where N is the number of patterns
%   types. Each cell contains the row, column, time active and pattern type
%   index for each occurence of the pattern type in that cell.
%
% PARAMS is a structure containing optional user settings:
%   params.minDurationSteps is the minimum duration (in time steps) of a 
%       pattern for it to be detected (default 1).
%   params.maxTimeGapSteps gives the maximum duration in time steps between
%       critical points (or synchrony/plane waves) for them to be counted
%       as the same pattern (default 0).
%   params.planeWaveThreshold gives the minimum order parameter for
%       activity to be considered a plane wave (must be between 0 and 1,
%       close to 1 indicates plane wave activity, default 1).
%   params.synchronyThreshold gives the maximum mean vector magnitude below
%       which all activity is considered synchony (default is two standard
%       deviations below the mean).
%   params.minEdgeDistance gives the minimum number of indices from the arrays
%       edge for critical points to be considered (default 0).
%   params.minCritRadius gives the minimum spatial radius for a critical
%       point to occupy for it to be counted, quantified by the winding
%       number and divergence/curl (default 1).
%   params.maxDisplacement gives the maximum number of grid spaces between
%       between critical points in consecutive patterns for them to be
%       counted as the same patterns (default 0.5).
%   params.combineNodeFocus is a boolean flag indicating that nodes and
%       foci should be treated as identical critical points, instead of
%       being split into sources/saddles and spirals (default false).
%   params.combineStableUnstable is a boolean flag indicating that stable
%       and unstable nodes and foci should be treated as identical critical
%       points types, instead of being split into sources and sinks
%       (default false).
% PARAMS can also be output, to keep track of any default values used.
%
%
% Rory Townsend, Aug 2018
% rory.townsend@sydney.edu.au


if any(size(vfx) ~= size(vfy))
    error('X and Y components of velocity field must be equal size!')
end

nt = size(vfx, 3);

[phi, v0, vdir] = orderParameter(vfx, vfy);
if exist('phase', 'var')
    rlength = nanmean(nanmean(exp(1i*phase)));
end
if ~exist('repNumber', 'var') || ~isnumeric(repNumber) || ~isscalar(repNumber)
    repNumber = [];
end
%% Check parameters and set default values if not supplied
if exist('params', 'var') && isstruct(params)
    inputParams = fieldnames(params);
else
    inputParams = [];
end

% Define list of default parameter values
checkParams = {'minDurationSteps', 'planeWaveThreshold', ...
    'synchronyThreshold', 'maxTimeGapSteps', 'minEdgeDistance', ...
    'minCritRadius', 'maxDisplacement', 'combineNodeFocus', ...
    'combineStableUnstable'};
defaultVals = [1 0.85 0.85 0 0 1 0.5 0 0];
isDefault = false(size(checkParams));

% Convert boolean values to doubles (temporarily)
if any(strcmp(inputParams, 'combineNodeFocus')) && ...
        islogical(params.combineNodeFocus)
    params.combineNodeFocus = double(params.combineNodeFocus);
end
if any(strcmp(inputParams, 'combineStableUnstable')) && ...
        islogical(params.combineStableUnstable)
    params.combineStableUnstable = double(params.combineStableUnstable);
end

% Set default values if any parameters are missing or invalid
for iparam = 1:length(checkParams)
    thisPar = checkParams{iparam};
    if ~any(strcmp(inputParams, thisPar)) || isempty(params.(thisPar)) ...
            || ~isnumeric(params.(thisPar)) || ~isscalar(params.(thisPar)) ...
            || params.(thisPar)<0
        params.(thisPar) = defaultVals(iparam);
        isDefault(iparam) = true;
    end
end

% Convert back to boolean
params.combineNodeFocus = params.combineNodeFocus == 1;
params.combineStableUnstable = params.combineStableUnstable == 1;

% Make sure that edge distance is at least as large as mininum radius
if params.minEdgeDistance < params.minCritRadius
    params.minEdgeDistance = params.minCritRadius;
end


%% Ask user for input on plane wave threshold if not supplied
% Work in progress, is not very helpful at this stage!
% figure
% if isDefault(strcmp(checkParams, 'planeWaveThreshold'))
%     if verLessThan('matlab', '9.2')
%         hist(phi, linspace(0.025, 0.975, 20))
%     else
%         histogram(phi, linspace(0,1,21))
%     end
%     xlabel('Order parameter')
%     ylabel('Counts')
%     title('Please specify plane wave threshold in command window')
%     disp('Specify order parameter minimum threshold for plane wave detection.')
%     disp('Value must be between 0 and 1, default 0.85.')
%     userThresh = input('Minimum order parameter: ');
%     if isempty(userThresh) || ~isnumeric(userThresh) || userThresh<0 || userThresh>0 
%         userThresh = 0.85;
%     end
%     params.planeWaveThreshold = userThresh;
% end
% 
%% Ask user for input on synchrony threshold if not supplied
% figure
% if isDefault(strcmp(checkParams, 'synchronyThreshold'))
%     if verLessThan('matlab', '9.2')
%         hist(rlength, 30)
%     else
%         histogram(rlength, 30)
%     end
%     axis tight
%     xlabel('Mean vector magnitude')
%     ylabel('Counts')
%     title('Please specify synchrony threshold in command window')
%     disp('Specify vector magnitude maximum threshold for synchrony detection.')
%     disp('Value must be between 0 and 1, default is 0.85.')
%     userThresh = input('Maximum vector magnitude: ');
%     if isempty(userThresh) || ~isnumeric(userThresh) || userThresh<0 || userThresh>0 
%         userThresh = 0;
%     end
%     params.synchronyThreshold = userThresh;
% end

%% Find plane waves and periods of synchrony
% Set up patterns array
pattTypes = {'Plane wave', 'Synchrony', 'Sink', 'Source', ...
    'Spiral-in', 'Spiral-out', 'Saddle'};
if params.combineNodeFocus
    if params.combineStableUnstable
        pattTypes = pattTypes([1:3, 7]);
        pattTypes{3} = 'Node/focus';
    else
        pattTypes = pattTypes([1:4, 7]);
    end
elseif params.combineStableUnstable
    pattTypes = pattTypes([1:3, 5, 7]);
    pattTypes{3} = 'Node';
    pattTypes{4} = 'Spiral';
end
allPatternLocs = cell(length(pattTypes), 1);

colNames = {'type', 'startTime', 'endTime', 'duration', 'startRow', ...
    'startCol', 'endRow', 'endCol', 'meanDisplacement'};
if ~isempty(repNumber)
    colNames = {colNames{:}, 'repetitionNumber'};
end
patterns = nan(nt, length(colNames));

% Plane waves
pwActive = phi>=params.planeWaveThreshold;
[pwStart, pwEnd, pwValid] = findRuns(pwActive, params.minDurationSteps, [], ...
    params.maxTimeGapSteps);
npw = length(pwStart);
patterns(1:npw, 1) = 1;
patterns(1:npw, 2:3) = [pwStart, pwEnd];
npat = npw;
% Find travel direction of plane waves
vdir = vdir(pwValid);
dirTimes = find(pwValid);
allPatternLocs{1} = cat(2, real(vdir), imag(vdir), dirTimes, ...
    ones(size(dirTimes)));

% Synchrony
% Only detect synchrony if phase data is input
if exist('phase', 'var')
    syActive = rlength(1:size(vfx,3)) >= params.synchronyThreshold;
    [syStart, syEnd] = findRuns(syActive, params.minDurationSteps, [], ...
        params.maxTimeGapSteps);
    nsy = length(syStart);
    patterns(npat + (1:nsy), 1) = 2;
    patterns(npat + (1:nsy), 2:3) = [syStart, syEnd];
    npat = npat + nsy;
    allPatternLocs{2} = zeros(0, 4);
end

%% Find all critical points
% Define map to quickly link critical point names to indices
if ~params.combineNodeFocus
    if ~params.combineStableUnstable
        critTypeInds = 1:5;
    else
        critTypeInds = [1 1 2 2 3];
    end
else
    if ~params.combineStableUnstable
        critTypeInds = [1 2 1 2 3];
    else
        critTypeInds = [1 1 1 1 2];
    end
end
ntypes = max(critTypeInds);
critTypeKeys = {'stableNode', 'unstableNode', 'stableFocus', ...
    'unstableFocus', 'saddle'};
ctype2ind = containers.Map(critTypeKeys, critTypeInds);

% Initialize cell array to store all valid critical points
critpointLocs = cell(nt, ntypes);
%%
for it = 1:nt
    % OPTIONAL: Skip time point if the system displays synchrony
    %if syActive(it)
    %    continue
    %end    
    % Find critical points (ignoring too close to edge)
    ivx = vfx(:,:,it);
    ivy = vfy(:,:,it);
    [rowcoords, colcoords, cptypes] = classifyCrit(ivx, ivy, ...
        params.minEdgeDistance);
    %%
    % Combine coordinates and add time
    allcoords = [rowcoords, colcoords, repmat(it, size(rowcoords))];
    cpIsValid = true(size(rowcoords));
    cpInds = zeros(size(cptypes));
 %%   
    % Loop over every critical point
    for icrit = 1:length(rowcoords)
        icoord = allcoords(icrit, :);
        cpInds(icrit) = ctype2ind(cptypes{icrit});
        
        % Calculate size of patterns using winding number
        if params.minCritRadius >= 1
            if strcmp(cptypes{icrit}, 'saddle')
                goodIndex = 1;
            else
                goodIndex = -1;
            end            
            % Iteratively calculate every winding number up to required
            % radius
            for irad = 1:params.minCritRadius
                index = windingNumberAngles(ivx, ivy, icoord, irad, 'both');
                if any(index ~= goodIndex)
                    cpIsValid(icrit) = false;
                    break
                end
            end                 
        end
    end
    
    % Save critical point types and locations
    cpInds = cpInds(cpIsValid);
    allcoords = allcoords(cpIsValid,:);
    
    for itype = 1:ntypes
        icoords = allcoords(cpInds==itype, :);
        critpointLocs{it, itype} = [icoords, (itype+2) * ...
            ones(size(icoords,1), 1)];
    end
end

%% Quick subfunction to find Euclidean distance between 2 coordinates
function d = findDist(pos1, pos2)
    d = sqrt(sum((pos1-pos2).^2, 2));
end

end
function [Flag,vehPos,indCUE,indDUE,indDUE2] = genCUEandDUE(d0, laneWidth, numLane, disBstoHwy, d_avg, numCUE,numDUE)
% First, generate traffic on the highway according to Spatial Poisson Process with
%      density determined by vehicle speed v (represented by avergage
%      inter-vehicle distance d_avg = 2.5*v).
% Then, generate DUE and CUE positions.

% Input: d0: highway length
%        d_avg: average inter-vehicle distance, d_avg = 2.5*v
%        Others should be evident from their names
% Output: Flag, 0/1, if the generated vehicles are not enough, flag = 1
%         vehPos, Nx2 matrix, storing vehicle coordinates
%         indCUE, numCUEx1, storing indices of CUEs in vehPos
%         indDUE, numDUEx1, storing indices of DUE transmitter
%         indDUE2, numDUEx1, storing indices of corresponding DUE receiver,
%                             closest to the DUEs stored in indDUE. 

% By Le Liang, Georgia Tech, Jan. 25, 2017

vehPos = []; % initilizer for vehicle position
indCUE = [];
indDUE = [];
indDUE2 = [];
Flag = 0;
%% generate all vehicle positions and store in vehPos
for ilane = 1:numLane
    npoints = poissrnd(2*d0/d_avg);
    pproc = (rand(npoints,1)*2-1)*d0; % horizontal coordinates
    pproc2 = [pproc, (disBstoHwy+ilane*laneWidth)*ones(length(pproc), 1)]; % [horizon vertical] coordinates
    vehPos = [vehPos; pproc2];
end
numVeh = size(vehPos, 1);
if numVeh < numCUE+2*numDUE
    Flag = 1; % the generated vehicles are not enough
    return; 
end
%% generate CUE and DUE positions
indPerm = randperm(numVeh);
indDUE = indPerm(1:numDUE); % randomly pick numDUE DUEs
indDUE2 = zeros(1,numDUE); % corresponding DUE receiver
for ii = 1 : numDUE
    % pair each element in indDUE with closet vehicle and store the
    % index in indDUE2
    minDist = 2*d0;
    tmpInd = 0;
    for iii = 1:numVeh
        if any(abs(iii-indDUE)<1e-6) || any(abs(iii-indDUE2)<1e-6) % iii in indDUE or indDUE2
            continue;
        end
        newDist = sqrt((vehPos(indDUE(ii), 1)-vehPos(iii,1))^2 + (vehPos(indDUE(ii), 2)-vehPos(iii,2))^2);
        if newDist < minDist
            tmpInd = iii;
            minDist = newDist;
        end
    end
    indDUE2(ii) = tmpInd; % the closest DUE pair
end

cntCUE = numDUE+1; % excluding those in indDUE
while cntCUE <= numVeh
    if any(abs(indPerm(cntCUE)-indDUE2)<1e-6) % element in indDUE2
        % do nothing
    else
        indCUE = [indCUE indPerm(cntCUE)];
    end
    cntCUE = cntCUE + 1;
    if length(indCUE) >= numCUE
        break
    end
end
indCUE = indCUE(:);
indDUE = indDUE(:);
indDUE2 = indDUE2(:);

end


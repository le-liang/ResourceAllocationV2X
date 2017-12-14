function [ assignment, minCapacity ] = maxMin( capacityMat )
%MAXMIN Algorithm for linear assignment problem which maximizes the minimum
%value of the assignment, rather than the max sum using Hungarian
%   Definition similar to [assignment, cost] = munkres(costMat)
% By Le Liang, Georgia Tech, July 29, 2016


[M, K] = size(capacityMat);
costMat1D = reshape(capacityMat, M*K, 1);
[sortVal] = sort(costMat1D, 'ascend');
minInd = 1;
maxInd = K*M;
assignment = ones(1,M);
while (maxInd - minInd) > 1
    mid = floor((minInd + maxInd)/2);
    tmpMat = capacityMat;
    for in = 1 : M
        for ik = 1 : K
            if tmpMat(in,ik) < sortVal(mid)
                tmpMat(in,ik) = 1;
            else tmpMat(in,ik) = 0;
            end
        end
    end
    [asgn, cost] = munkres(tmpMat);
    if cost > 0 
        maxInd = mid;
    else
        minInd = mid;
        assignment = asgn;
    end
end
minCapacity = sortVal(minInd);

end


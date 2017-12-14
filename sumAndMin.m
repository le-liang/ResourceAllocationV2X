function [ sumVal, minVal ] = sumAndMin( capMat, assignment )
%SUMANDMIN returns the sum and min of the chosen elements assignment in
% the matrix capMat
% Input: capMat, the computed capacity matrix from main function
%        assignment, a vector containing the column indices of each row 
%           corresponding to the chosen elements
% Output: sumVal, the sum value of all chosen elements
%         minVal, the minimum value of all chosen elements
%  By Le Liang, Georgia Tech, July 29, 2016

vec = zeros(length(assignment), 1);
M = size(capMat, 1);

for irow = 1 : M
    vec(irow) = capMat(irow, assignment(irow));
end

sumVal = sum(vec);
minVal = min(vec);


end


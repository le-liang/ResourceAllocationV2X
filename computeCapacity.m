function [ output_args ] = computeCapacity( a, b )
% compute the CUE capacity according to the closed form expression
% a/((a-b)*log2)*(exp(1/a)*expint(1/a) - exp(1/b)*expint(1/b)) according to
% Lemma 2 in the paper

% By Le Liang, Georgia Tech, Jan. 25, 2017

if a>=(1/700) && b>=(1/700)
    output_args = a/((a-b)*log(2))*(exp(1/a)*expint(1/a) - exp(1/b)*expint(1/b));
elseif a<(1/700) && b<(1/700)
    output_args = a/((a-b)*log(2))*(a - b);
elseif b < (1/700)
    output_args = a/((a-b)*log(2))*(exp(1/a)*expint(1/a) - b); % exp(x)*expint(x)= 1/x when x is big, interesting finding
elseif a < (1/700)
    output_args = a/((a-b)*log(2))*(a - exp(1/b)*expint(1/b));
end

end


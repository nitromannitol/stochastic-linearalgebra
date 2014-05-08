function [transformMat] = LHC(n,d,epsilon)
%implements JLT based off of Ping Li, Trevor Hastie, Kenneth Church, 2006
%
%TO DO: Merge with ACH and introduce S variable for sparsity
%
% NOTE: WHAT ARE MARGINAL NORMS?
%
%Syntax: transformMat = LHC(n,d, epsilon)
%
%
% Inputs:
%   n = number of points
%   d = dimension of each point
%	epsilon = distortion constant, smaller = slower but more accurate
%
%
%
% Outputs:
%    transformMat = the projection matrix
%
%
%------------------------------------------------------------------

if(d > n)
    n
    d
	error('Must have more points that dimension')
end

%gurantee pairwise distance preservation
%with probability n^{-b} 
beta = 0.1;
k = ceil((4 + 2*beta)*log(n) / (epsilon^2/2 - epsilon^3/3) );
if( k > d)
    error('Adjust epsilon or beta, k is too large')
end

%Create the very sparse R matrix as described in the paper 
S = sqrt(d);
probability = (1/(2*S));

R = zeros(k,d);
decision = rand(k,d);
R(decision < probability) = sqrt(3);
R(decision > (1- probability)) = -sqrt(3);
R = R/sqrt(k);

transformMat = R;

end
function [transformMat] = ACH(n,d,epsilon)
%implements JLT based off of Achlioptas, 2001
%
%
%
%Syntax: transformMat = ACH(n,d, epsilon)
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

%Create the R matrix as described in the paper 
R = zeros(k,d);
decision = rand(k,d);
R(decision < 1/6) = sqrt(3);
R(decision > 5/6) = -sqrt(3);
R = R/sqrt(k);

transformMat = R;

end
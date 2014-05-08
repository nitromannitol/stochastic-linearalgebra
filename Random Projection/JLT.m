function [transformMat] = JLT(n,d,epsilon)
%implements a standard JLT 
%
%
%Syntax: transformMat = JLT(n,d,epsilon)
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
if(1/sqrt(epsilon) > d)
    1/sqrt(epsilon)
    d
	error('1/sqrt(epsilon) must be less than input dimension')
end

%Set target dimension
k = ceil(4*log(n)/(epsilon^2 - epsilon^3/3));
if(k > d)
    error('Epsilon too low')
end

%transformMat = zeros(k,d);
%decision = rand(k,d);
%transformMat(decision < 1/2) = 1;

%G = randn(d,k);
%[Q,R] = qr(G);
%transformMat = Q';

tmp = rand(d,1); tmp = 2*(tmp>1/2)-1;
tmp = diag(tmp);
F = exp(2*pi*i*([0:d-1]'*[0:d-1]/d));
aux = ceil(rand(k,1)*k);
F = F(aux,:);
aux
size(F)
size(tmp)
transformMat = F*tmp;





end
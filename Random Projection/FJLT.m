function [transformMat] = FJLT(n,d,epsilon,p)
%implements FJLT based off of Ailon, Chazello Approximate Nearest Neighbors and FJLT
%
% The complexity for p =1 is O(dlog d + epsilon^{-3} log^2 n)
%
%Syntax: transformMat = FJLT(A,B, k)
%
%
% Inputs:
%   n = number of points
%   d = dimension of each point
%	epsilon = distortion constant, smaller = slower but more accurate
%	p = the dimensional of the functional space, i.e., L_p	    
%
%
%
% Outputs:
%    transforMat = the projection matrix
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

%Set some large enough constant, 5 for now
%test out different constants later 
c = 5;

%Set target dimension
k = ceil(c*log(n)/epsilon^2);
if(k > d)
    error('Epsilon too low')
end


%Create the P, H, and D matricies as described in section 2 on page 2 of
%the paper
%
%
%P is a k-by-d matrix whose elements are independent mixtures of 0 with an unbiased normal disbirution 
%of variance 1/q, where q = min(theta(epsilon^(p-2) log^p(n)/d, 1))
% P_ij - is normal(0,inv(q)) with probabiilty q and - with probability 1 - q
%
%first define q
q = min( (log2(n)^p * epsilon^(p-2)/d),1);


%With probability q, P_ij is normal(0,q^-1), 0 otherwise
P = rand(k,d);
%With probabilty 1-q, P_ij is 0
P = (P < q);
%With probability q, P_ij is normal(0,q^-1)
P = (1/q)*randn(k,d) .* P;



%H is a d-by-d normalized Hadamard matrix,
%H_ij = d^{-1/2} (-1)^{<i-1,j-1>}
H = 1/d^(0.5)*hadamard(d);



%D is a d-by-d diagonal matrix, where each D_ii is drawn independently from {-1,1} with 
%probability 1/2
diagonalEntries = ones(d,1);
diagonalEntries(rand(d,1) > 0.5) = -1;
D = diag(diagonalEntries);

transformMat = P*H*D;

end
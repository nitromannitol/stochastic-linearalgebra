function [kLargestElems] = getKLargestMatrixMultiply(A, B, k, projection_type)
%This returns the k-largest elements of the matrix product A*B
%
% This is pretty slow and has been proved to have a large
% variance, i.e., not accurate for most cases
% see: http://cseweb.ucsd.edu/~akmenon/HonoursThesis.pdf
% section 5.1
%
%Syntax: kLargestElements = getKLargestMatrixMultiply(A,B, k)
%
%
% Inputs:
%   A,B = two matrices
%   k = the number of largest in magnitude elements you would like to get from A*B
%   projectionType = type of projection you want to use, it defaults 
%   to FJLT for now
%   These are the options 
%   'None' - Uses a zero matrix, used to check baseline error rates
%   'FJLT' - uses fast johnson lindenstrauss transform 
%   'JLT' - uses standard johnson lindenstrauss transform
%   'ACH'  - uses the sparse projection by achlioptas, 2001, with s = 3
%   'LHC' - uses the very sparse projection by Li,Hastie,Church (2006)
%  
%   Note: in terms of speed of computing the projection matrix, FJLT and JLT are 
%   the best by a factor of 8 or so
%
%   
%  
%
% Outputs:
%    kLargestElems = an array containing the kLargest elements
%
%------------------------------------------------------------------


if(nargin < 4)
    projection_type = 'FJLT';
end

if(strcmp(projection_type,'FJLT'))
    dim = size(A,2);
     %Pad A and B to the nearest power of 2
     nextPower = 2^nextpow2(dim);
     A = [A zeros(size(A,1),nextPower - dim)];
     B = [B; zeros(nextPower - dim, size(B,2))];    
end

n = size(A,1) + size(B,2);

if(k > n)
    k
    n
    error('k must be less than the number of data points')
end


d = size(A,2);

if(d ~= size(B,1))
    error('Matrix multiplication mismatch')
end


epsilon = 0.4;
%
%TO DO: Figure out what p does
%
p = 1; 


%Use the maximum likelihood estimate of 
%the inner product approximating the distribution
%by a normal distribution and estimating
% the inner product using the MLE, which is
% just the average 


% Number of times to repeat
delta = 0.1;
%numRepeats = ceil(1/delta); 
numRepeats = 1;

%Keep track of total
C_total = zeros(size(A,1),size(B,2));


for i = 1:numRepeats

    %Choose projection
    if(strcmp(projection_type,'FJLT'))
        project_mat = FJLT(n,d,epsilon,p);
    elseif(strcmp(projection_type,'JLT'))
        project_mat = JLT(n,d,epsilon);
    elseif(strcmp(projection_type,'ACH'))
        project_mat = ACH(n,d,epsilon);
    elseif(strcmp(projection_type,'LHC'))
        project_mat = LHC(n,d,epsilon);
    elseif(strcmp(projection_type,'None'))
        project_mat = zeros(k,d);
    end
     
    %k
    %size(project_mat)
    %norm(project_mat*project_mat' - eye(size(project_mat,1)))
    eig(project_mat*project_mat')
    
    %
    %Project A and B and calculate inner product in lower dimensional space
    C_proj = project(A,B,project_mat);
    
    %Keep track of the total
    C_total = C_total + C_proj;
end

%Take the average
C_total = C_total/numRepeats;

%Order the elements in the lower dimensional space
[sortedVals, sortedIndex] = sort(C_total(:), 'descend');

k_proj_index = sortedIndex(1:k);
[I,J] = ind2sub(size(C_total),k_proj_index);


%Compute the k largest elems in the higher-dimensional space
%This takes O(k*n^2) time
kLargestElems = zeros(k,1);
for curr_elem = 1:k
    i = I(curr_elem);
    j = J(curr_elem);
    kLargestElems(curr_elem) = A(i,:)*B(:,j);
end







end
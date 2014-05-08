%test script

%Test speed and manually test accuracy
clear all


%initialize 2 random matricies with mean 0 and std deviation 1000
A = 1000.*randn(2000,1000);
B = 1000.*randn(1000,2000);


%Test getKLargestMatrixMultiply
k = 10; 
tic
kLargest = getKLargestMatrixMultiply(A,B,k, 'JLT')
toc
kLargest = sort(kLargest(:), 'descend')

%actually compute the k largest elements with brute force
tic
C = A*B;

[sortedValues] = sort(C(:),'descend');                                               
maxVals = sortedValues(1:k)
toc

%%
tic 
d = 1000;
for i = 1:10
     project_mat = ACH(n,d,epsilon);
end
toc






%%
%Estimate the error rate

clear all; 
totalError = 0;
numIters = 50;

for i = 1:numIters
    
%pick 2 random matricies again 
A = 50.*randn(500,256);
B = 50.*randn(256,500);



%Test getKLargestMatrixMultiply
k = 10; 
kLargest = getKLargestMatrixMultiply(A,B,k, 'FJLT');
kLargest = sort(kLargest,'descend');


%actually compute the k largest elements with brute force
C = A*B;

[sortedValues] = sort(C(:),'descend');                                               
maxVals = sortedValues(1:k);


%Just some error metric rate to test
totalError =  totalError + norm(maxVals - kLargest)/norm(maxVals);
end


%note that a benchmark error is ~1, done with just picking the first k
%elems from C
totalError = totalError/numIters


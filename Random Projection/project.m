function [C_proj] = project(A,B,projectMat)
%Returns the projected inner product of elements in
% A and B using the projectedMat
%
%
%Syntax: [A_proj, B_proj] = project(A,B,projectMat)
%
% Inputs:
%    A = matrix of whose rows you want to project 
%    B = a matrix whose columns you want to project
%    projectMat = a projection matrix,e.g., a JT
%
%  Outputs:
%   C_proj = a matrix of the projected rows in the low_dim space
% 
%
%------------------------------------------------------------------

%Project A and B down to lower dim
A_proj = projectMat*A';
B_proj = projectMat*B;

%Computer inner product in lower dimensional space
C_proj = A_proj'*B_proj; 



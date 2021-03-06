function [buffer] = putVectorInBuffer(in,buffer, n)
%putVectorInBuffer Puts a whole vector into a buffer from position N to
%N+length(in)

% IN is a column
% N is the first position that will be filled, from left to righ 

% BUFF = [0 0 0 0 0]';
% N = 3;
% IN = PAST<[ 1 2 ]'>PRESENT
% PUTIN = [0 0 1 2 0];

len = length(buffer);
indexC = mod(n-1,len)+1;
buffer(indexC: indexC+length(in)-1,  1) = in;

end


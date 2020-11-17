% Increased efficiency by new input
function [buffer] = putVectorInBufferV3(in,buffer, bufferLength, n)
%putVectorInBuffer Puts a whole vector into a buffer from position N to
%N+length(in)

% IN is a column
% N is the first position that will be filled, from left to righ 

% BUFF = [0 0 0 0 0]';
% N = 3;
% IN = PRESENT <[ 1 2 ]'>PAST
% PUTIN = [1 2 0 0 0];

indexC = mod(n-1,bufferLength)+1;
indexL = indexC+length(in)-1; % Last position
buffer(indexC: indexL,  1) = in;



if indexL > bufferLength
    
    a = buffer(bufferLength+1:indexL);
        buffer( 1: (indexL  -   bufferLength) ) = a;
    buffer = buffer(1:bufferLength);



end


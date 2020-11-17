% Increased efficiency by new input
function [buffer] = putVectorInBufferV4(in,buffer, BufferLength, n)
% Increased efficiency by not changing the size of the buffer under no
% condition (as did in V3)
%putVectorInBuffer Puts a whole vector into a buffer from position N to
%N+length(in)

% IN is a column
% N is the first position that will be filled, from left to righ 

% BUFF = [0 0 0 0 0]';
% N = 1;
% IN = PRESENT <[ 1 2 ]'>PAST
% PUTIN = [1 2 0 0 0];

indexC = mod(n-1,BufferLength)+1;
indexL = indexC+length(in)-1; % Last position

    if indexL < BufferLength

        buffer(indexC: indexL,  1) = in;
        
    else

        buffer(indexC: BufferLength,  1) = in(1:BufferLength-indexC+1);
        buffer(1: abs(indexL - BufferLength),  1) = in((BufferLength-indexC+2): end);

    end

end


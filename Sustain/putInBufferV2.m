% Increased efficiency by new input
function [buffer] = putInBufferV2(in, buffer, bufferLength, n)
%CIRCULARBUFFER Summary of this function goes here
%   Detailed explanation goes here
% in SHOULD BE A COLUMN
indexC = mod(n-1,bufferLength)+1;
buffer(indexC,1) = in;

end


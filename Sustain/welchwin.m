function w = welchwin(n)
% WELCHWIN  Welch (parabolic) window
%
% Usage:
%   w = welchwin(n) returns the n-point parabolic window
%         n: scalar window length
%
%         w: column vector of window values
% 
% Example:
%   t = (0:9)';
%   s = rand(1,10)';
%   n = length(s);
%   w = welchwin(n);
%   s2 = s.*w;
%   plot(t, s, 'ko', t, s2, 'r-')
%   legend('original data', 'windowed data')
%
% See also: RECTWIN

% v0.1 (Nov 2012) by Andrew Davis (addavis@gmail.com)
% modelled after RECTWIN (built-in)
% some window definitions: http://paulbourke.net/miscellaneous/windows/

narginchk(1,1)
assert(isscalar(n), 'n must be a scalar')

j = (0:n-1)';
w = 1 - (2*j/n - 1).^2;    % parabolic window

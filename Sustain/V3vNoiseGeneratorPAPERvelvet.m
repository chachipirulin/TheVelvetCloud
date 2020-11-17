
function [vn, se, k, NumberOfImpulses ] = V3vNoiseGeneratorPAPERvelvet(Ls,Fs, Nd, DecayConstant)

%VNOISEGENERATOR Version v2
% Consist on ramdomly spaced unitary bipolar impulses [-1,0,1]

% EXAMPLE --> Density = 0.1; Samples = 44100;

vn = zeros(Ls, 1); % <PREALLOC>

% Grid size
 Td = Fs / Nd; %Nd: The pulse density, or the average number of nonzero impulses per secon

% How many impulses in this VN?
NumberOfImpulses = floor(Ls / Fs * Nd ); 
m = (1:NumberOfImpulses)';

% Positions at where impulses are located?
k = floor( m * Td + rand(NumberOfImpulses, 1) * (Td - 1));

    % make sure our locations are less than NoiseSamples size
k = k(k < length(vn));    

% Sign (-1,1) for each impulse
sign = 2 * round(rand(NumberOfImpulses,1)) - 1; 

% We need to ramdomize the starting position of the first impulse (and conseuently, the others)
% k = k + round(rand * Ls * Density * 10);

% Assing the sign to those positions
vn(k) = sign(1:length(k));

    %% Computation of Se
% r3(m) is a random gain between 0.5 and 2.0
r3 = rand(NumberOfImpulses,1) * 1.5 + 0.5;

se = exp(-DecayConstant .* m) .* sign .* r3;

     % Make sure we are not taking more
se = se(1:length(k));

% Make sure #Samples out = #Samp in
% vn = vn(1:Ls);

end


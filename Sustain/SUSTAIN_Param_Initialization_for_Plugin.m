% Input - Definitions - Preinitializations -  Random preallocations

%__________________ Input _______________________%
[in,Fs] = audioread('chords.wav');
PostStart = 13231; % So we start FROM this position to acess past 
in = in(PostStart:end,1); 
in_mono = in(:,1); % mono-1, 1]
N = length(in_mono);
ButtonTakeNewSnippet = true; % When you press the pedal, take a nex snippet

%__________________ Definitions ___________________%
SustainThreshold = 0.6; %[0 1] gain has to be higher than so sustain acts
NumberOfVoicesAverage = 32;  % Select average number of voices per sample OR grain size: 30 / second
Nd = 80; % Impulses per second in Velvet noise, THE HIGHER, THE NOISER!!!!! THE LOWER, THE "BEATING"!
DecayConstant = 0.001; % Decay constant for the impulses in vn
WindowsLength = round(3 * 10e-2 * Fs); % 30ms, dangelo example last graph

BufferLength = WindowsLength;
Po = 2e-5; % Pressure reference level, NOT USED IN V3

%__________________ Preinitializations ___________________%

out_mono = zeros(length(in_mono), 1);
inSnippetBuffer = zeros(WindowsLength, 1); % Input buffer to hold a snippet
outSnippet = zeros(WindowsLength, 1); % Not really used, just to hold values for plotting
outSamples = zeros(N, 1);
playbackVoice = zeros(WindowsLength, 1);
LsZerosVector = zeros(WindowsLength, 1);
Vn = zeros(N, 1);
SnippetSample = zeros(N, 1);
WindowsStatic = zeros(WindowsLength, 1);

%__________________ Random preallocations for positions and signs ___________________%

signsArray = [     1     1    -1     1    -1    -1     1    -1    -1     1     1     1     1     1     1     1     1    -1    -1     1    -1    -1     1    -1    -1    -1    -1     1]';
positionsArray = [87         148        1157        2140        2345        2633        3235        3245        3410        3481        3538        3963        4346        4362 4507        4786        4858        5337        6082        7755        7980        9020        9060        9137        9552       10200       10324       12356]';
% already sorted
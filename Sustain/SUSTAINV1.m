% Comments at bottom

classdef SUSTAINV1 < audioPlugin         

%% PROPERTIES #############################################################
    properties
        Fs = 44100;
        One_time_trigger = 0;
        PostStart = 13231; % So we start FROM this position to acess past 
        init = false;       % (Re)Initialise variables
    end
% Initialized when plugin is created
    properties (Access = private)

%__________________ Definitions ___________________%
SustainThreshold = 0.5; %[0 1] gain has to be higher than so sustain acts
NumberOfVoicesAverage = 32;  % Select average number of voices per sample OR grain size: 30 / second
Nd = 80; % Impulses per second in Velvet noise, THE HIGHER, THE NOISER!!!!! THE LOWER, THE "BEATING"!
DecayConstant = 0.001; % Decay constant for the impulses in vn
WindowsLength = 13230; %round(3 * 10e-2 * Fs); % 30ms, dangelo example last graph
%ButtonTakeNewSnippet = true; % When you press the pedal, take a next snippet--> THIS IS SUBTITUTED BY IF BUFFER IS FULL OR NOT, IF YOOU PRESS, IT IS "empty"
IsSnippetBufferFull = false; % We will may need many IN instances to fill the 30ms snippet
SamplesAlreadyInBuffer = 0;
BufferLength = 13230;
Po = 2e-5; % Pressure reference level, NOT USED IN V3

%__________________ Preinitializations ___________________%

% out_mono = zeros(length(in_mono), 1);
inSnippetBuffer = zeros(13230, 1); % Input buffer to hold a snippet
outSnippet = zeros(13230, 1); % Not really used, just to hold values for plotting
playbackVoice = zeros(13230, 1);
LsZerosVector = zeros(13230, 1);
Vn = zeros(1, 1);
WindowsStatic = zeros(13230, 1);

%__________________ Random preallocations for positions and signs ___________________%

% signsArray = [      -1     1    -1    -1     1    -1    -1     1     1     1     1     1     1     1     1    -1    -1     1    -1    -1     1    -1    -1    -1    -1     1]';
% positionsArray = [      148        1157        2140        2345        2633        3235        3245        3410        3481        3538        3963        4346        4362 4507        4786        4858        5337        6082        7755        7980        9020        9060        9137        9552       10200       10324       12356]';
% % already sorted
    end
      
    methods
%% MAIN LOOP #############################################################
function out = process(p,in)      
    
                % mono normalized
            in_mono = in(:,1);% ./ max(abs(in(:,1)) + 10e-17); % So we don't divide by 0 when silent = NaN
            
                % Initialise the out-vector
            out_mono = zeros(size(in_mono));        
      
        % If we have a snippet to convolve, we proceed
    if p.IsSnippetBufferFull == true

        signsArray = [  1    -1     1    -1    -1     1    -1    -1     1     1     1     1     1     1     1     1    -1    -1     1    -1    -1     1    -1    -1    -1    -1     1]';
        positionsArray = [      148        1157        2140        2345        2633        3235        3245        3410        3481        3538        3963        4346        4362 4507        4786        4858        5337        6082        7755        7980        9020        9060        9137        9552       10200       10324       12356]';
% already sorted

            %% Update Equation
        for n = 1 : length(in_mono) % 1 y 2
            
            % We already have our snippet buffer filleD

                % Velvet noise for this time
            p.Vn = getVnoiseSample(p.Fs, p.Nd);

                if p.Vn  ~= 0 % then it is an impulse, update buffers
                  signsArray = [ p.Vn; signsArray ];
                  positionsArray = [ 1 ; positionsArray ]; % Now is 1, Take the past(est) out
                end

            % We avance one, one delay for the postision of impulses
        positionsArray = positionsArray + 1;
        
            if positionsArray(end) > p.WindowsLength    
                    % If an impulse position has reached the end      
                positionsArray = [positionsArray(1:end - 1)];
                         % bye bye last
                 signsArray = [signsArray(1:end - 1)]; 
                 
            end

            % Multiply the impulses at positions by the incoming Input snippet
        p.playbackVoice = sum( p.inSnippetBuffer(positionsArray) .* signsArray);

            % Susbtraction
        out_mono(n,1) =  p.playbackVoice;% Mono, I know
        
        
        end
%         max(in(:,1))
    out = [out_mono zeros(length(in_mono), 1)]; % Here you are your stereo
    
        % If a value of our snippet was higher than thres, reinit the
        % buffer
            if (max(in_mono)) > p.SustainThreshold
                p.IsSnippetBufferFull = false;
                display("We HAVE REACHED THE THRESLHOLD")
            end
    else
            % Buffer is not full, IN=OUT
        out = in; % TODO: meanwhile buffer is loading, we should output also a convlution, not IN=OUT!
        
        % Fill the buffer meanwhile

            % Put this little snippet in buffer
        p.inSnippetBuffer = putVectorInBufferV2(in_mono, p.inSnippetBuffer, p.BufferLength, length(in_mono));
        
            % Update how full is the buffer
        p.SamplesAlreadyInBuffer = length(in_mono) + p.SamplesAlreadyInBuffer;

                % If buffer is full
            if p.SamplesAlreadyInBuffer > p.WindowsLength

                    % Then we can window our definite snippet to
                    % convolve
                p.inSnippetBuffer = p.inSnippetBuffer .* p.WindowsStatic;
                    % And buffer is full
                p.IsSnippetBufferFull = true;
                p.SamplesAlreadyInBuffer = 0;

            end
        end
end
        
%% METHODS ################################################################

                % Buffer functions
            function [out] = accessBufferIndexes(buffer, delayPosition, n)
                len = length(buffer);
                indexD = mod(n - delayPosition - 1,len) + 1;
                out = buffer(indexD,1);
            end
            function [buffer] = putVectorInBufferV2(in,buffer, bufferLength, n)
            %putVectorInBuffer Puts a whole vector into a buffer from position N to
            %N+length(in)

            % IN is a column
            % N is the first position that will be filled, from left to righ 

            % BUFF = [0 0 0 0 0]';
            % N = 3;
            % IN = PAST<[ 1 2 ]'>PRESENT
            % PUTIN = [0 0 1 2 0];

            indexC = mod(n-1,bufferLength)+1;
            buffer(indexC: indexC+length(in)-1,  1) = in;

            end

                % VN generator: Generates one sample of Velvet noiseat density Nd per second (Fs)
            function [VnoiseSample] = getVnoiseSample(Fs, Nd)


            VnoiseSample = 0;

                if rand < Nd/Fs
                    VnoiseSample = round(rand) * 2 - 1;
                end

            end

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
            end

                 % Run at the begining and no more
            function reset (p)

                if  p.One_time_trigger == 0 % Create the WELCH WINDOW one time! 

                        p.WindowsStatic = welchwin(p.WindowsLength);
                        p.One_time_trigger = ++p.One_time_trigger;
                        p.IsSnippetBufferFull = false;

                end

            end

        end
        
    end




%% PLOT RESULTS ################################################################


%% RUN PLUGIN ################################################################
% validateAudioPlugin SUSTAINV1.m
% audioTestBench SUSTAINV1.m
% generateAudioPlugin SUSTAINV1

%% Properties TRICKS: https://se.mathworks.com/help/audio/ug/tips-and-tricks-for-plugin-authoring.html
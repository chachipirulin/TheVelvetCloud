% Comments at bottom
 
    classdef SUSTAINV7 < audioPlugin         
 
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
SustainThreshold = 0.6; %[0 1] gain has to be higher than so sustain acts
NumberOfVoicesAverage = 32;  % Select average number of voices per sample OR grain size: 30 / second
Nd = 300; % Impulses per second in Velvet noise, THE HIGHER, THE NOISER!!!!! THE LOWER, THE "BEATING"!
DecayConstant = 0.001; % Decay constant for the impulses in vn
WindowsLength = 13230; %round(3 * 10e-2 * Fs); % 30ms, dangelo example last graph
%ButtonTakeNewSnippet = true; % When you press the pedal, take a next snippet--> THIS IS SUBTITUTED BY IF BUFFER IS FULL OR NOT, IF YOOU PRESS, IT IS "empty"
IsSnippetBufferFull = false; % We will may need many IN instances to fill the 30ms snippet
SamplesAlreadyInBuffer = 0;
BufferLength = 13230;
Po = 2e-5; % Pressure reference level, NOT USED IN V3
FirstRise = 0;
howFullIsBuffer = 0;
ReachingMax = 0;
playbackVoiceGain = 0.1;
sustainKnobGain = 0.2;
% maximunsArray = zeros(1024*,1); Y% 18ms --> 795 samples
 
%__________________ Low Pass coefs ffor convolving snippet _________________%
    % fcut
Fcut= 5000;
Wn =  5000 / 44100;   
order = 3;
    % Coeffs of butterworth design
b = zeros(1, 3);
a = zeros(1, 3);
 
%__________________ Preinitializations ___________________%
 
% out_mono = zeros(length(in_mono), 1);
inSnippetBuffer = zeros(13230, 1); % Input buffer to hold a snippet
CopyOf_inSnippetBuffer = zeros(13230, 1);
outSnippet = zeros(13230, 1); % Not really used, just to hold values for plotting
playbackVoice = zeros(1, 1);
LsZerosVector = zeros(13230, 1);
Vn = zeros(1, 1);
VnBuffer = zeros(1, 13230);
WindowsStatic = zeros(13230, 1);
playbackVoiceSumBuff = zeros(13230, 1);
arrayOfImpulses = zeros(1024, 1);
 
%__________________ Random preallocations for positions and signs ___________________%
 
signsArray = [  1  -1  -1   1  -1  1    -1    -1     1    -1     1     1     1     1     1        1      -1    -1     1    -1    -1     1    -1    -1    -1    -1     1  zeros(1,100) ]';
positionsArray = [  148     500   750  1000 1157        2140             2633        3235        3481        3538        3963        4346        4362 4507        4786        4858        5337        6082        7755        7980        9020        9060        9137        9552       10200       10324       12356  zeros(1,100)]';
% already sorted and zero padded
TotalImpulses = 27;
 
    end
      
    methods
%% MAIN LOOP #############################################################
function out = process(p,in)      
    
    %% BEFORE RUNNING, RUN THIS LINE BELOW
    % p.ReachedAMax = 0;
                % mono normalized
            in_mono = in(:,1);%./ max(abs(in(:,1)) + 10e-17); % So we don't divide by 0 when silent = NaN
            
                % Initialise the out-vector
            out_mono = zeros(size(in_mono));        
            
                % Initialise the NORMALIZED out-vector
            %out_mono_norm = zeros(size(in_mono));  % NOT USED      
      
        % If buffer is NOT full or gain > threshold...
            if p.IsSnippetBufferFull == false % you want to "convolve" with this incomming snippet
 
                    % We put the chunk in the buffer
               p.inSnippetBuffer = putVectorInBufferV4(in_mono, p.inSnippetBuffer, p.BufferLength, p.SamplesAlreadyInBuffer + 1 );
 
                    % Update how full is buffer
               p.SamplesAlreadyInBuffer = p.SamplesAlreadyInBuffer + length(in_mono);  
                    
                    % If buffer is completed with the chunks,
                if p.SamplesAlreadyInBuffer > p.BufferLength
                
                        % It is full,
                    p.IsSnippetBufferFull = true; 
                    disp("bufferfull");
                        % windowize it,
                    p.inSnippetBuffer = p.inSnippetBuffer .* p.WindowsStatic;
                    
                        % Lowpass filter it
                     p.inSnippetBuffer = filter(p.b, p.a, p.inSnippetBuffer);    
                        
                        % Duplicate, actually, we will use the copy for convolution
                     p.CopyOf_inSnippetBuffer = p.inSnippetBuffer .* p.playbackVoiceGain; % InSnippetBuff will be as high as our ThresholdLow
                        
                        % and reinit
                     p.SamplesAlreadyInBuffer = 0;   
                end
            end
 
            %% Update Equation
        for n = 1 : length(in_mono) % 1 y 2
            
 
                % Velvet noise for this time
            p.Vn = getVnoiseSample(p.Fs, p.Nd);
                
                if p.Vn  ~= 0 % then it is an impulse, update buffers
                    
                    % DUMMY buffer to ckeck everything is right
%                  p.VnBuffer = putVectorInBufferV3(p.Vn , p.VnBuffer, p.BufferLength, n );
                    
                    % Update total impulses
                 p.TotalImpulses = p.TotalImpulses + 1;
 
                    % Put the Vn at this present position
                 p.signsArray = [p.Vn ; p.signsArray(1 : end-1) ];
                 p.positionsArray = [1; p.positionsArray(1 : end-1) ]; % Now is 1
                 
                end
                
                
            % We avance one, one delay for the postision of impulses
        p.positionsArray = p.positionsArray + 1;
        
            
            if p.positionsArray(p.TotalImpulses) > p.WindowsLength
                      % If an impulse position has reached the end   
                 p.signsArray(p.TotalImpulses) = 0; % set its sign to zero
                 p.positionsArray(p.TotalImpulses) = 0; % set its position to zero
                 
                    % Update total impulses
                p.TotalImpulses = p.TotalImpulses - 1;
                 
            end
            
            %DUMMY to check the #numberof impulses we have at each sample
%          p.arrayOfImpulses(n) = p.TotalImpulses;
%          figure;title("#ofImpulses per sample");plot(p.arrayOfImpulses); stem(mean(p.arrayOfImpulses>0));
         
            % Multiply the impulses by the COPY of the Input snippet
        p.playbackVoice = sum( p.CopyOf_inSnippetBuffer(p.positionsArray(1:p.TotalImpulses)) .* p.signsArray(1:p.TotalImpulses)); 
 
        out_mono(n,1) =  p.playbackVoice(1, 1); % TODO: When we also substract sounded better, so we can check...
 
        end
 
            %Plot the #of Impulsesin this input chunk and its mean
        
            % Mix Sustain & Input
        out = p.sustainKnobGain .* out_mono + (1 - p.sustainKnobGain) .* in_mono;
       
        out = [out zeros(length(in_mono), 1)]; % Here you are your stereo
    
        % If a value of our snippet was higher than thres, add 1 to  first
        % rise
    if any((in_mono) > p.SustainThreshold )
        p.ReachingMax = 1; % DUMMY to check
        p.FirstRise = p.FirstRise + 1;
        
                % If actually it was a first rise (=1), then take the
                % snippet
            if p.FirstRise == 1
                p.IsSnippetBufferFull = false;
            end
            
%         disp("We HAVE REACHED THE THRESLHOLD: " + p.FirstRise);
 
    elseif any((in_mono) < p.SustainThreshold )
        
        p.ReachingMax = 0; % DUMMY to check
        p.FirstRise = 0;
%         disp("We FIRST rise, reinit buffer: " + p.FirstRise);
    end
%             
            

        close all;
%     figure; hold on; 
%     subplot(331); stem(p.IsSnippetBufferFull); title("p.IsSnippetBufferFull");
%     subplot(332); plot(p.CopyOf_inSnippetBuffer); hold on; title("p.CopyOf_inSnippetBuffer");
%     subplot(336); plot(out) ; title("out(n) = in* (1-knob)+ sus*knob"); % If you plot when it's reseted, YOU WILL GET 0!
%     subplot(334);plot(in_mono); title("in_mono");
%     subplot(335);plot(p.inSnippetBuffer); title(" p.inSnippetBuffer");
%     subplot(337);stem(p.ReachingMax); title("p.ReachingMax");
%     subplot(338);stem(p.FirstRise); title("p.FirstRise");
%      subplot(333);plot(p.sustainKnobGain .* out_mono ); title("p.sustainKnobGain .* out_sustain ");
%     subplot(339);stem(p.CopyOf_inSnippetBuffer(p.positionsArray(1:p.TotalImpulses)) .* p.signsArray(1:p.TotalImpulses)); title("inSnippet(positions)*signs(n)");
%  
 
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
 
                        p.SustainThreshold = 0.8; %[0 1] gain has to be higher than so sustain acts
                        p.WindowsStatic = welchwin(p.WindowsLength);
                        p.inSnippetBuffer = zeros(13230, 1);
                        p.CopyOf_inSnippetBuffer = zeros(13230, 1);
                        p.One_time_trigger = ++p.One_time_trigger;
                        p.IsSnippetBufferFull = false;
                        p.SamplesAlreadyInBuffer = 0;
                        p.playbackVoice = 0;
                        p.Nd = 120;
                        p.VnBuffer = zeros(1, 13230);
p.signsArray = [  1  -1  -1   1  -1  1    -1    -1     1    -1     1     1     1     1     1        1      -1    -1     1    -1    -1     1    -1    -1    -1    -1     1  zeros(1,100) ]';
p.positionsArray = [  148     500   750  1000 1157        2140             2633        3235        3481        3538        3963        4346        4362 4507        4786        4858        5337        6082        7755        7980        9020        9060        9137        9552       10200       10324       12356  zeros(1,100)]';
% already sorted and zero padded
                        p.TotalImpulses = 27;
                        p.playbackVoiceGain = 0.2;
                        p.sustainKnobGain = 0.5;

                        %  _____ Filter definitionfor convolving snippet ________
                       [p.b,p.a] = butter(p.order, p.Wn, 'low');
                end
 

            end
 
        end
        
    end
 
 
%% PLOT RESULTS ################################################################
% 
%     close all;
%     figure; hold on; 
%     subplot(331); stem(p.IsSnippetBufferFull); title("p.IsSnippetBufferFull");
%     subplot(332); plot(p.CopyOf_inSnippetBuffer); hold on; title("p.CopyOf_inSnippetBuffer");
%     subplot(333); plot(p.playbackVoice); title("playbackVoice"); % If you plot when it's reseted, YOU WILL GET 0!
%     subplot(334);stem(in_mono); title("in_mono");
%     subplot(335);plot(p.inSnippetBuffer); title(" p.inSnippetBuffer");
%     subplot(337);stem(p.ReachingMax); title("p.ReachingMax");
%     subplot(338);stem(p.FirstRise); title("p.FirstRise");
%     subplot(336);plot(out(:,1)); title("out");
 
 
%% RUN PLUGIN ################################################################
% validateAudioPlugin SUSTAINV7.m
% audioTestBench SUSTAINV7.m
% generateAudioPlugin SUSTAINV7
 
%% VERSION COMMENTS
 
%   - Convolving Snippet is triggered just in the first rise (copy from V3)
%   - LP filter for convolving snippet to avoid Hfreq
%   - Listening to the mean of voices:
%       Nd40 : mean= ARTIFACTS
%       Nd80 : mean=
%       Nd120 : mean=
%       Nd150 : mean=
%       Nd200 : mean=
%
 
 
%% Properties TRICKS: https://se.mathworks.com/help/audio/ug/tips-and-tricks-for-plugin-authoring.html

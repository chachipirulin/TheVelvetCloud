% Comments at bottom

% QUESTIONS/ PROBLEM:

%     - MegaBuffer changes size at some point and it is not supossed to do
%       so! when generateAudioPlugin FDN4V3
%     - Output is 0 ...

classdef FDN4V4 < audioPlugin         

%% PROPERTIES #############################################################
    properties
        Fs = 44100;
        Trigger_Only_One_Time = 0;
        init = false;       % (Re)Initialise variables
    end
% Initialized when plugin is created
    properties (Access = private)


%__________________ Delay Lengths _________________________________________

   % Max delay of 999 ms
maxDelay = ceil(.999 * 44100);
d1 = fix(.0297*44100); d2 = fix(.0371*44100); d3 = fix(.0411*44100);d4 = fix(.0437*44100);

%___________________________ Buffers ______________________________________
BufferSize = 2000; % d4
% Buffers must be able to hold the highest delay ("d4")
   % Initialize all buffers
buffer1 = zeros(1900,1); 
buffer2 = zeros(2000,1); 
buffer3 = zeros(2000,1); 
buffer4 = zeros(2000,1);

MegaBufferSize = 9000;
megaBuffer = zeros(9000,1); 

%__________________ Fractional delay param ___________________
% LFO parameters
rate1 = 0; amp1 = 0; rate2 = 0; amp2 = 0; rate3 = 0; amp3 = 0; rate4 = 0; amp4 = 0; 
    
%_________________  Feedback matrices ___________________

% A = 0.5 * [ -1  1 1 1 1 1 1 1;...   % Scattering matrix: % HOUSEHOLD
%             1  -1  1 1 1 1 1 1;...
%             1 1  -1  1 1 1 1 1;...
%             1 1 1  -1  1 1 1 1;...
%             1 1 1 1  -1  1 1 1;...
%             1 1 1 1 1  -1  1 1;...
%             1 1 1 1 1 1  -1  1;...
%             1 1 1 1 1 1 1   -1];

a11   = -1 ; a12 =   1 ; a13 =   1 ; a14 =   1 ; 
a21   =  1 ; a22 =   -1; a23 =   0 ; a24 =   1 ; 
a31   =  1 ; a32 =   1 ; a33 =  -1 ; a34 =   1 ; 
a41   =  1 ; a42 =   1 ; a43 =  -1 ; a44 =  -1 ;

%___________________ Output/Input coefficients ___________________________

NoiseSamples = 8200;
PulseDensity = 1000; % pulses per second (Fs)
NoiseDuration = 10; % IN MILISECONDS!! paper says = LS!!!!
Density = 1000;
DensityPercentage = 0.1;
DecayConstant = 0.1;

g =.4; % Output delay lines Gain, to control reverb time USED DUMMY INSTEAD OF B-C

% ___________________ Feedback holding variables ___________________
fb1 = 0; fb2 = 0; fb3 = 0; fb4 = 0;
outB1 = 0; outB2 = 0; outB3 = 0; outB4 = 0; 
InputLengthHolder = 0;
%__________________ Noise buffers - IN RESET FUNCTION Computed only one time!!! ___________________

vn1,  NumberOfImpulses; vn2, vn3, vn4;
k1 = zeros(185, 1);
se1 = zeros(185, 1);
k2 = zeros(185, 1);
se2 = zeros(185, 1);
k3 = zeros(185, 1);
se3 = zeros(185, 1);
k4 = zeros(185, 1);
se4 = zeros(185, 1);
idx_process = 0; % To count how many times the process is runned in the test bench
One_time_trigger = 0; 
Cur_idx =0; % Cyrrent index

    end
    
    
    methods
%% MAIN LOOP #############################################################
function out = process(p,in)      
    
                % mono 
            in_mono = in(:,1);% ./ max(abs(in(:,1)) + 10e-17); % So we don't divide by 0 when silent = NaN
            
                % Initialise the stereo out-vector
            out_mono = zeros(size(in_mono));      
            
                % Fill megaBuffer putVectorInBufferV2(in,buffer, bufferLength, n)
            p.megaBuffer = putVectorInBufferV4(in_mono, p.megaBuffer, p.MegaBufferSize, length(in_mono) + p.InputLengthHolder); % Write locations          

            % Past will be contained in the MegaBuffer, initialized with
            % silence. Its past will be accessed at K positions
            
            %% Update Equation
        for n = 1 : length(in_mono) % 1 y 2

                % Reinit INPUT coefficients
            p.outB1 = 0; p.outB2 = 0; p.outB3 = 0; p.outB4 = 0; 

         % x * Se
                    % outB1
                x_se1 = sum(accessBufferIndexes(p.megaBuffer, p.k1, p.InputLengthHolder + n).* p.se1);
                p.outB1 = x_se1;
                    % outB2
                x_se2 = sum(accessBufferIndexes(p.megaBuffer, p.k2, p.InputLengthHolder + n).* p.se2);
                p.outB2 = x_se2; 
                    % outB3
                x_se3 = sum(accessBufferIndexes(p.megaBuffer, p.k3, p.InputLengthHolder + n).* p.se3);
                p.outB3 = x_se3;
                    % outB4
                x_se4 = sum(accessBufferIndexes(p.megaBuffer, p.k4, p.InputLengthHolder + n).* p.se4);
                p.outB4 = x_se4;
                
                
         % inDL: Sum outB1 with fb (feedback coeffs) for each delay line
         
                % Combine input with feedback for respective delay lines
            inDL1 = p.outB1 + p.fb1; inDL2 = p.outB2 + p.fb2; inDL3 = p.outB3 + p.fb3; inDL4 = p.outB4 + p.fb4;

         % outDL: Parallel output delay lines
            [outDL1, p.buffer1] = modDelay(inDL1,p.buffer1, p.Fs, n + p.idx_process * length(in_mono) , p.d1, p.amp1, p.rate1);
            [outDL2, p.buffer2] = modDelay(inDL2,p.buffer2,p.Fs, n + p.idx_process * length(in_mono) , p.d2,p.amp2,p.rate2);
            [outDL3, p.buffer3] = modDelay(inDL3,p.buffer3,p.Fs, n + p.idx_process * length(in_mono) , p.d3,p.amp3,p.rate3);
            [outDL4, p.buffer4] = modDelay(inDL4,p.buffer4,p.Fs, n + p.idx_process * length(in_mono) , p.d4,p.amp4,p.rate4);

          % fb: Calculate feedback coefficients entering the matrix (including crossover)
            p.fb1 = p.g*(p.a11*outDL1 + p.a21*outDL2 + p.a31*outDL3 + p.a41*outDL4);
            p.fb2 = p.g*(p.a12*outDL1 + p.a22*outDL2 + p.a32*outDL3 + p.a42*outDL4);
            p.fb3 = p.g*(p.a13*outDL1 + p.a23*outDL2 + p.a33*outDL3 + p.a43*outDL4);
            p.fb4 = p.g*(p.a14*outDL1 + p.a24*outDL2 + p.a34*outDL3 + p.a44*outDL4);

                % Combine parallel paths
            out_mono(n,1) = 0.1*(outDL1 + outDL2 + outDL3 + outDL4); % Mono, I know
%             
%                 close all
% figure;  hold on;
% subplot(331); plot(in_mono); title("in_mono");
% subplot(332);plot(out_mono); title("Out");
% subplot(333);plot(p.megaBuffer); title("megaBuffer");
% subplot(334);plot( accessBufferIndexes(p.megaBuffer, p.k1, n + length(in_mono))); title("accessBufferIndexes(p.megaBuffer, p.k1...))");
% subplot(335);stem(inDL1, 'r*');title("inDL1");
% subplot(336);plot( accessBufferIndexes(p.megaBuffer, p.k3, n + length(in_mono))); title("accessBufferIndexes(p.megaBuffer, p.k3...))");
% subplot(337);plot(p.buffer1); title("p.buffer1");
% subplot(338);stem(outDL1, 'r*'); title("outDL1");
% subplot(339);stem(x_se2, 'r*'); title("x_se3");

            % SET BREAKPOINT HERE AND LOOK HOW BUFFER 1 GETS FILLLED!!!
        end
        
    out = [out_mono zeros(length(in_mono), 1)]; % Here you are your stereo
    
        % Update the lentgh position holder of the input
    p.InputLengthHolder = p.InputLengthHolder + length(in_mono);
%         p.InputLengthHolder = mod(p.InputLengthHolder + length(in_mono),
%         length(p.megaBuffer)); % silvin

                % Update current index 
        p.idx_process = p.idx_process + 1;
                
%     close all
% figure;  hold on;
% subplot(331); plot(in_mono); title("in_mono");
% subplot(332);plot(out_mono); title("Out");
% subplot(333);plot(p.megaBuffer); title("megaBuffer");
% subplot(334);plot( accessBufferIndexes(p.megaBuffer, p.k1, n + length(in_mono))); title("accessBufferIndexes(p.megaBuffer, p.k1...))");
% subplot(335);stem(inDL1, 'r*');title("inDL1");
% subplot(336);plot( accessBufferIndexes(p.megaBuffer, p.k3, n + length(in_mono))); title("accessBufferIndexes(p.megaBuffer, p.k3...))");
% subplot(337);plot(p.buffer1); title("p.buffer1");
% subplot(338);stem(outDL1, 'r*'); title("outDL1");
% subplot(339);stem(x_se2, 'r*'); title("x_se3");

end
        
%% METHODS ################################################################

            % Buffer functions
        function [out] = accessBufferIndexes(buffer, delayPosition, n)
            len = length(buffer);
            indexD = mod(n - delayPosition - 1,len) + 1;
            out = buffer(indexD,1);
        end
        function [out,buffer] = circularBuffer(in,buffer, delay, n)


        len = length(buffer);
        indexC = mod(n-1,len)+1;

        indexD = mod(n-delay - 1,len) + 1;
        out = buffer(indexD,1);

        buffer(indexC,1) = in;
        end
        
            % VN generator
                %  inputs: NoiseSamples, Fs, PulseDensity, DecayConstant); 
        function [vn, se, k, NumberOfImpulses ] = V3vNoiseGeneratorPAPERvelvet(Ls,Fs, Nd, DecayConstant)

            vn = zeros(Ls, 1); % <PREALLOC>

            % Grid size
             Td = Fs / Nd; %Nd: The pulse density, or the average number of nonzero impulses per secon

            % How many impulses in this VN?
            NumberOfImpulses = floor(Ls / Fs * Nd ); 
            m = (1:NumberOfImpulses)';

            % Sign (-1,1) for each impulse
            sign = 2 * round(rand(NumberOfImpulses,1)) - 1; 

            % Positions at where impulses are located?
            k = round( m * Td + rand(NumberOfImpulses, 1) * (Td - 1));

            % We need to ramdomize the starting position of the first impulse (and conseuently, the others)
            % k = k + round(rand * Ls * Density * 10);

            % Assing the sign to those positions
            vn(k) = sign;

            % Computation of Se
            % r3(m) is a random gain between 0.5 and 2.0
            r3 = rand(NumberOfImpulses,1) * 1.5 + 0.5;

            se = exp(-DecayConstant .* m) .* sign .* r3;

        end

            % ModDelay
        function [out,buffer] = modDelay(in,buffer,Fs,n, delay,depth,rate)

                   % Calculate time in seconds for the current sample
                t = (n-1) / Fs;
                fracDelay = depth * sin(2*pi*rate*t);
                intDelay = floor(fracDelay);
                frac = fracDelay - intDelay;

                   % Determine indexes for circular buffer
                len = length(buffer);
                indexC = mod(n - 1,len) + 1;   % Current index
                indexD = mod(n - delay - 1 + intDelay,len) + 1; % Delay index
                indexF = mod(n - delay - 1 + intDelay+1,len) + 1; % Fractional index


                % Delay index indexF = mod(n–delay–1+intDelay+1,len) + 1; % Fractional index
                out = (1 - frac) * buffer(indexD,1) + (frac) * buffer(indexF,1);
                   % Store the current output in appropriate index
                buffer(indexC,1) = in;
                end

             % Run at the begining and no more
        function reset (p)
            
            if  p.One_time_trigger == 0 % Create the noises one time! 
                    [p.vn1, p.se1, p.k1, p.NumberOfImpulses] = V3vNoiseGeneratorPAPERvelvet(p.NoiseSamples, p.Fs, p.PulseDensity, p.DecayConstant); 
                    [p.vn2, p.se2, p.k2] = V3vNoiseGeneratorPAPERvelvet(p.NoiseSamples, p.Fs, p.PulseDensity, p.DecayConstant); 
                    [p.vn3, p.se3, p.k3] = V3vNoiseGeneratorPAPERvelvet(p.NoiseSamples, p.Fs, p.PulseDensity, p.DecayConstant); 
                    [p.vn4, p.se4, p.k4] = V3vNoiseGeneratorPAPERvelvet(p.NoiseSamples, p.Fs, p.PulseDensity, p.DecayConstant);

                    p.One_time_trigger=++p.One_time_trigger;

            end
            
%             p.InputLengthHolder = 0;
            p.megaBuffer = zeros(9000,1); 
            p.InputLengthHolder = 0;

            p.buffer1 = zeros(2000,1); % TODO:NEED TO RESIZE THIS TO MAX DELAY LINES
            p.buffer2 = zeros(2000,1); 
            p.buffer3 = zeros(2000,1); 
            p.buffer4 = zeros(2000,1); 
              p.idx_process = 0;
        end
        
                    % Not used
        function p = FDN4
            
            % Call the initialise function on creation of the plugin
            % UNUSED!!!!!!
            initialise (p);
        
        end
        
    end

end

%% PLOT RESULTS ################################################################
% close all
% figure;  hold on;
% subplot(331); plot(in_mono); title("In");
% subplot(332);plot(out_mono); title("Out");
% subplot(333);plot(p.megaBuffer); title("megaBuffer");
% subplot(334);plot( accessBufferIndexes(p.megaBuffer, p.k4, n + length(in_mono))); title("accessBufferIndexes(p.megaBuffer, p.k4...))");
% subplot(335);plot(x_se2); title("x_se2");
% subplot(336);plot(x_se3); title("x_se3");

%% RUN PLUGIN ################################################################
% validateAudioPlugin FDN4V4.m
% audioTestBench FDN4V4.m
% generateAudioPlugin FDN4V4

%% VERSION COMMENT ################################################################

%   - STATE: WORKING
 % TODO:NEED TO RESIZE THIS TO MAX DELAY LINES
%   - Listen to output: silvin.wav

%% Properties TRICKS: https://se.mathworks.com/help/audio/ug/tips-and-tricks-for-plugin-authoring.html


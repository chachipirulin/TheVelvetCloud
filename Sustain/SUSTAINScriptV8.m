% Comments at bottom

% ___________________________ Parameters __________________________%

SUSTAIN_Param_Initialization_for_Plugin;

% _________________ Initialization (from functions) _______________%

    % Window to get snippet
WindowsStatic = welchwin(WindowsLength);

% ________________________ Sample loop _____________________________%

    for n = 1:N
           
        % CONVOLVING SNIPPET:
            if ButtonTakeNewSnippet == true % you want to "convolve" with this incomming snippet
                
                inSnippetBuffer = putVectorInBufferV2(in_mono((PostStart-BufferLength+1):PostStart), inSnippetBuffer, BufferLength, n);
                inSnippetBuffer = inSnippetBuffer .* WindowsStatic;
                
                    % until you press again the button
                ButtonTakeNewSnippet = false; 
                
            end

                % Velvet noise for this time
            Vn(n) = getVnoiseSample(Fs, Nd);

                if Vn(n)  ~= 0 % then it is an impulse, update buffers
                 signsArray = [ Vn(n); signsArray ];
                 positionsArray = [ 1 ; positionsArray ]; % Now is 1, Take the past(est) out
                end

            % We avance one, one delay for the postision of impulses
        positionsArray = positionsArray + 1;
        
            if positionsArray(end) > WindowsLength    
                    % If an impulse position has reached the end      
                positionsArray = [positionsArray(1:end - 1)];
                         % bye bye last
                 signsArray = [signsArray(1:end - 1)]; 
                 
            end

            % Multiply the impulses by the incoming Input snippet
        playbackVoice(n) = sum( inSnippetBuffer(positionsArray) .* signsArray);

            % Susbtraction
        outSamples(n,1) =  playbackVoice(n);% If you substract, same sound
        
    end

       % Normalization
    outSamples = outSamples / max(abs(outSamples));

%% NOT USED FOR PLUGIN ####################################################

% ________________________ Sound results _____________________________%

       % Write & listen audio
%     audiowrite('chordsDANGELOsustainV8.wav', outSamples, Fs );     
%     soundsc(outSamples + in_mono .* 0.5, Fs )

% ________________________ Plot results _____________________________%

%     close all
%     figure; 
%     subplot(331); plot(inSnippetBuffer); hold on; title("inSnippetBuffer");
%     subplot(332); plot(SnippetSample); hold on; title("SnippetSample(n)");
%     subplot(333); plot(playbackVoice); hold on; title("playbackVoice"); % If you plot when it's reseted, YOU WILL GET 0!
%     subplot(334);plot(outSnippet); title("outSnippet");
%     subplot(334);plot(Vn); title("Vn(n)");
%     subplot(335);plot(outSamples); title("outSamples");
%     subplot(336);plot(); title("nothing");
%     subplot(336);plot(outSamples); title("outSamples");
% 
% soundsc(outSamples, Fs)


%% DBs CONVERSION #############  OUT OF ORDER  ############################
% Comment from Stefano: 
% "Plus, you don't actually need the compressor part if you ask me, but
% still you can have it if you want." --> then it should work without this
% (?)
  
% [LinearGain]  = snippetToVolume(in(4000:5000), playbackVoice(4000:5000), Po);

%% PLOTS FOR DBs ##########  OUT OF ORDER  #################################

%     % Graphs
% close all; figure; hold on; 
% subplot(331); plot(in_mono); title("Input_D_C")
% subplot(332); plot(out_mono); title("Output")
% subplot(333); plot(w); title("Window")
% subplot(334); plot(inWindowed); title("Input_W_i_n_d_o_w_e_d")
% 
% subplot(335); plot(inSnippetBuffer); title("voiceBuffer(n)")
% subplot(336); plot(dBdiff); title("dBdiff")
% subplot(337); plot(out_substracted); title("Output substracted")
% subplot(338); plot(dBplaybackVoice); title("dB playbackVoice")
% subplot(339); plot(linearDiff); title("linearDiff")

%% _________________________________________________________________________
% @author: Carmen ML
% @version comment: 
%         - "Copy" of DANGELOsustainScriptV8
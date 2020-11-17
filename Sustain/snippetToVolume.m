function [LinearDiff] = snippetToVolume(inSnippet,playbackVoiceSnippet, Po)
% Comment from Stefano: 
% "Plus, you don't actually need the compressor part if you ask me, but
% still you can have it if you want." --> then it should work without this



        % Full wave rectification
    in = abs(inSnippet);     
    playbackVoice = abs(playbackVoiceSnippet);

        % Avoid negative log!
     in = in + 1e-17;
     playbackVoice = playbackVoice + 1e-17;  
        
        % to dB, Max-min linear = [0 1] ----- Max-min dBs [-140  93.9794] (aprox with Po)
    dBin = sum(20 * log10(in/Po))/ numel(in);
    dBplaybackVoice = sum(20 * log10(playbackVoice / Po)) / numel(playbackVoice);
    
        % Difference
    dBdiff = dBin - dBplaybackVoice ;
    
        % Back to linear
     LinearDiff = 20.^((dBdiff + 10*log10(Po) )/ 10);

        % Plot results
%          close all; figure; hold on; 
% subplot(321); plot(in); title("in")
% subplot(322); plot(playbackVoice); title("playbackVoice")
% subplot(323); stem(dBin); title("dBin")
% subplot(324); stem(dBplaybackVoice); title("dBplaybackVoice")
% subplot(325); stem(abs(dBdiff)); title("dBdiff")
% subplot(326); stem(abs(LinearDiff)); title("LinearDiff")


end


classdef ModSigGen
    % ModSigGen is a static class for generate modulation signal.
    % signal time is [0, param.T].
    enumeration
        CW, LFM, PSK, FSK, FD
    end
    
    methods(Static)
        function sig = generate(param) 
            t = linspace(0, param.T, param.T * param.FS);
            switch param.type
                case ModSigGen.CW 
                    sig = ModSigGen.cw(param.amp, param.freq, param.phase, t);
                case ModSigGen.LFM
                    sig = ModSigGen.lfm(param.amp, param.freq, param.phase, t, param.rate);
                case ModSigGen.PSK
                    sig = ModSigGen.psk(param.amp, param.freq, param.phase, t, param.code, param.encodingNum);
                case ModSigGen.FSK
                    sig = ModSigGen.fsk(param.amp, param.freqList, param.phase, t);
                case ModSigGen.FD
                    sig = ModSigGen.fd(param.amp, param.freqList, param.phase, t);
                otherwise
                    sig = zeros(size(t));
            end
        end
        
        function y = cw(am, f, ph, t)
           y = am .* cos(2 * pi * f * t + ph); 
        end
        
        function y = lfm(am, f, ph, t, k)
            y = am .* cos(2 * pi * (f * t + k * t.^2 ./ 2) + ph);
        end
        
        function y = psk(am, f, ph, t, code, encodingNum)
            n = length(code);
            encodingPhase = [0, pi];
            if encodingNum == 4
                encodingPhase = [0, pi / 2, pi, pi * 3 / 2];
            end
            deltaDecodeWidth = t(end) / n; 
            y = zeros(size(t));
            for i = 1 : n
                rect = ones(size(t)) .* (abs((t - (i - 0.5) * deltaDecodeWidth) ./ deltaDecodeWidth) <= 0.5);
                y = y + am * rect .* cos(2 * pi * f * t + encodingPhase(code(i) + 1) + ph);
            end
        end
        
        function y = fsk(am, flist, ph, t)
            n = length(flist);
            deltaDecodeWidth = t(end) / n; 
            y = zeros(size(t));
            for i = 1 : n
                rect = ones(size(t)) .* (abs((t - (i - 0.5) * deltaDecodeWidth) ./ deltaDecodeWidth) <= 0.5);
                y = y + am * rect .* cos(2 * pi * flist(i) * t + ph);
            end
        end
        
        function y = fd(am, flist, ph, t)
            n = length(flist);
            y = zeros(size(t));
            for i = 1 : n
                y = y + am * cos(2 * pi * flist(i) * t + ph);
            end
        end
    end
end


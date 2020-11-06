function name = getTypeName(type)
    switch type
        case ModSigGen.CW 
            name = "CW";
        case ModSigGen.LFM
            name = "LFM";
        case ModSigGen.PSK
            name = "PSK";
        case ModSigGen.FSK
            name = "FSK";
        case ModSigGen.FD
            name = "FD";
        otherwise
            name = "NONE";
    end
end


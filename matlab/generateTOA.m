function toa = generateTOA(params, toa_range, pri_type)
% generate toa by pri_type, which between toa_range

    float_rate = 0.01; reject_rate = 0.1;
    
    % chosse pri generator type
    switch pri_type
        case PRIType.Stable
            f = @generateStablePRI;
        case PRIType.Wobble
            f = @generateWoblePRI;
        case PRIType.Sin
            f = @generateSinPRI;
        case PRIType.Stagger
            f = @generateStaggerPRI;
        case PRIType.Slip
            f = @generateSlipPRI;
        otherwise
            f = @generateStablePRI;
    end
    
    % generate toa
    tBegin = toa_range(1); tEnd = toa_range(2);
    toa = [];
    ind = 1;
    toa(ind) = tBegin;
    while toa(end) <= tEnd
        pri = f(params);
        toa(ind + 1) = toa(end)  + (1 + (2 * rand - 1) * float_rate) * pri;
        ind = ind + 1;
    end
    toa(end) = [];
    nTOA = length(toa);
    
    % reject a part of toa
    rejectIndList = randperm(nTOA, round(reject_rate * nTOA));
    toa(rejectIndList) = [];
    
    % pri generate function
    function pri = generateStablePRI(params)
        pri = params(1);
    end

    function pri = generateWoblePRI(params)
        pri = params(1);
        wobble_rate = params(2);
        pri = pri * (1 + sqrt(1 / 3) * randn * wobble_rate);
    end

    function pri = generateStaggerPRI(params)
        k = mod(length(toa) - 1, length(params));
        pri = params(k + 1);
    end

    function pri = generateSlipPRI(params)
        if length(toa) < 2
            pri = params(1);
        else
            pri = toa(end) - toa(end - 1) + params(2);
            if pri > params(2)
                pri = params(1);
            end
        end
    end

    function pri = generateSinPRI(params)
        mid_val = params(1);
        am = params(2);
        cycle = params(3);
        ph = params(4);
        pri = mid_val + am * mid_val * sin(2 * pi * length(toa) / cycle + ph);
    end
end


        

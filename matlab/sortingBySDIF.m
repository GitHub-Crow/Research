function signal = sortingBySDIF(toa, pri_range)
    %% set parameters
    MAX_PULSE_LOSS_RATIO = 0.5;
    DECAY = 0.8;
    BIN_COEFF = 0.6;
    minPRI = pri_range(1); maxPRI = pri_range(2);
    TOLERANCE = 0.02 * (maxPRI - minPRI);
    LEVEL = 5; % maximum level 
    MIN_EVENTS = 5; % minimim events count need for a squence search 

    %% CDIF
    diff = 1;
    % record(i, :) is [pri, begin, end] represent toa(end) - toa(begin) =
    % toa. signal(i) is {pri, seq}, seq is result sequence and pri is
    % primary pri.
    record = [];
    signal = struct('pri', {}, 'seq', {});
    TOA_END = toa(end);
    
    while diff < LEVEL
        record = [];
        buildKLevelHist(diff, toa);
        record = sortrows(record);
        % build hist by slip windows
        Lind = 1; Rind = 1; ind = 1;
        nRecord = length(record);
        priRecord = []; % record prime pri in the window
        indRecord = struct('begin', {}, 'end', {}); % record beign and end index of the window

        % detect pri and count valid pulse
        while ind <= nRecord
            while Rind <= nRecord && record(ind, 1) + TOLERANCE >= record(Rind, 1) % move right side of window
                Rind = Rind + 1;
            end
            while record(ind, 1) - TOLERANCE > record(Lind, 1)  % move left side of window
                Lind = Lind + 1;
            end
            % valid duration is [Lind, Rind)
            priRecord(end + 1) = record(ind, 1);
            indRecord(end + 1).begin  = Lind;
            indRecord(end).end = Rind - 1;
            ind = ind + 1;
        end
    
        nPRI = length(priRecord);
        th = MAX_PULSE_LOSS_RATIO * nPRI * exp(-priRecord ./ (BIN_COEFF * priRecord(end)));
        actualCount = zeros(1, nPRI);
        for k = 1 : nPRI
            actualCount(k) = indRecord(k).end - indRecord(k).begin + 1;
        end
        % debug info image
        close all; figure ;
        plot(priRecord, actualCount, 'LineWidth', 1.8); hold on
        plot(priRecord, th, 'LineWidth', 2.2); hold on
        ylabel('count', 'FontSize', 12, 'FontName', 'cambria');
        xlabel('PRI / us', 'FontSize', 12, 'FontName', 'cambria');
        legend('actual count', 'threshold', ...
               'FontSize', 12, 'FontName', 'cambria');
        % judge and search sequence
        isFindSeq = false;
        
        % only check peak PRI
        potentialPRIInd = find(actualCount >= th);
        [~, peakPRIInd] = findpeaks(actualCount);
        potentialPRIInd = intersect(potentialPRIInd, peakPRIInd);
        
        % search by actual event count 
        combinedArr = sortrows([actualCount(potentialPRIInd)', potentialPRIInd'], 'descend');
        potentialPRIInd = combinedArr(:, 2);
        potentialPRIInd = potentialPRIInd';
        nPInd = length(potentialPRIInd);
        % if peack not exceed threshold then do subharmonic checking
        [~, maxCountPRIInd] = max(actualCount);
        peakPRI = priRecord(maxCountPRIInd);
        if nPInd && ~ismember(maxCountPRIInd, potentialPRIInd)
            subharmonicTolerance = ceil(potentialPRIInd(1) / peakPRI) * TOLERANCE;
            if mod(priRecord(potentialPRIInd(1)), peakPRI) <= subharmonicTolerance ||...     % example : mod(15, 7)
               peakPRI - mod(priRecord(potentialPRIInd(1)), peakPRI) <= subharmonicTolerance % example : mod(21, 7)
                isFindSeq = searchSeq(maxCountPRIInd);
            end
        else
            if diff ~= 1 || nPInd == 1
                for pInd = potentialPRIInd
                    isFindSeq = searchSeq(pInd);
                    if isFindSeq, break ; end
                end
            end
        end
        if isFindSeq
            diff = 1;
            record = [];
            MAX_PULSE_LOSS_RATIO = MAX_PULSE_LOSS_RATIO * DECAY;
        else
            diff = diff + 1;
        end
    end

    

    function isSuccessful = searchSeq(k)
       pri = priRecord(k);
       indArr = [];
       seqTolerance = 0.4 * TOLERANCE;
       for j = indRecord(k).begin : indRecord(k).end
           ind_1st = record(j, 2); % first event index
           ind_2st = record(j, 3); % second event index
           toa_seq_begin = toa(ind_1st);
           tmpIndArr = [ind_1st, ind_2st];
           eventK = 2;
           while true
               idealToa =  toa_seq_begin + eventK * pri;
               if idealToa > TOA_END, break ; end
               [actualToa, ind_] = searchCloseTarget(toa, idealToa);
               if abs(actualToa - idealToa) <= eventK * seqTolerance && ...
                  ind_ ~= tmpIndArr(end)
                   tmpIndArr(end + 1) = ind_;
                   eventK = eventK + 1;
               else
                   break ;
               end
           end
           % judge whether search sequence successfully
           if eventK > MIN_EVENTS
               indArr = [indArr, tmpIndArr];
           end
       end
       indArr = unique(indArr);
       indArr = sort(indArr);
       if length(indArr) < MIN_EVENTS
           isSuccessful = false;
           return ;
       end
       isSuccessful = true;
       signal(end + 1).pri = pri;
       signal(end).seq = toa(indArr);
       toa(indArr) = [];
    end

    function buildKLevelHist(k, inArg)
        N = length(inArg);
        for i = 1 : N - k
            pri = inArg(i + k) - inArg(i);
            record(end + 1, :) = [pri, i, i + k];
        end
    end

    function [localPRI, localK] = searchCloseTarget(data, inArg)
        L = 1; R = length(data);
        while R - L > 1
            M = ceil((L + R) / 2);
            if data(M) < inArg
                L = M;
            elseif data(M) > inArg
                R = M;
            else
                localPRI = data(M);
                localK = M;
                return ;
            end
        end
        if inArg - data(L) < data(R) - inArg
            localPRI = data(L);
            localK = L;
        else
            localPRI = data(R);
            localK = R;
        end
    end
end
function signal = sortingByCDIF(toa, pri_range)
    %% set parameters
    % one level and two level threshold coefficient 
    % it is clear that wrong sorting extract more toa than actual
    % situation, which will pose bad influence on subsequent sorting, so
    % set DECAY to decrease threshold for correct sorting futher. 
    C1 = 0.7; C2 = 0.4; DECAY = 0.8;
    minPRI = pri_range(1); maxPRI = pri_range(2);
    TOLERANCE = 0.01 * (maxPRI - minPRI);
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
        
        % debug info image
        close all; figure ;
        actualCount = zeros(1, length(priRecord));
        for k = 1 : length(priRecord)
            actualCount(k) = indRecord(k).end - indRecord(k).begin + 1;
        end
        plot(priRecord, actualCount, 'LineWidth', 1.8); hold on
        plot(priRecord, C1 * TOA_END ./ priRecord, 'LineWidth', 2.2); hold on
        plot(priRecord, C2 * TOA_END ./ priRecord, 'LineWidth', 2.2); hold on
        ylabel('count', 'FontSize', 12, 'FontName', 'cambria');
        xlabel('PRI / us', 'FontSize', 12, 'FontName', 'cambria');
        xlim([minPRI, inf]);
        legend('actual count', 'first threshold', 'second threshold', ...
               'FontSize', 12, 'FontName', 'cambria');
           
        % only search peak PRI
        [~, peakIndArr] = findpeaks(actualCount);
        isFindSeq = false;
        for k = peakIndArr
            onePRI = priRecord(k);
            if onePRI < minPRI, continue ; end
            if onePRI > maxPRI, break ; end
            if TOA_END / MIN_EVENTS < onePRI, break ; end % there is no more events
            
            % find pri which has maximum actual events count and whose
            % value meet double of onePRI
            doubleEventCnt = max(actualCount(abs(priRecord - 2 * onePRI) <= TOLERANCE));
            oneEventCnt = actualCount(k);
            if oneEventCnt > C1 * TOA_END / onePRI && ...
               doubleEventCnt > C2 * TOA_END / (2 * onePRI)
                isFindSeq = searchSeq(k);
                break ;
            end
            
        end
        if isFindSeq
            diff = 1;
            record = [];
            C1 = C1 * DECAY;
            C2 = C2 * DECAY;
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
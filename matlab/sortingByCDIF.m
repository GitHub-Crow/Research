function signal = sortingByCDIF(toa, pri_range)
    %% set parameters
    % one level and two level threshold coefficient 
    % it is clear that wrong sorting extract more toa than actual
    % situation, which will pose bad influence on subsequent sorting, so
    % set DECAY to decrease threshold for correct sorting futher. 
    C1 = 0.8; C2 = 0.6; DECAY = 0.6;
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
        x = []; y = [];
        for k = 1 : length(priRecord)
            if priRecord(k) < minPRI, continue ; end
            x(end + 1) = priRecord(k);
            y(end + 1) = indRecord(k).end - indRecord(k).begin + 1;
        end
        plot(x, y, 'LineWidth', 1.8); hold on
        plot(x, C1 * TOA_END ./ x, 'LineWidth', 2.2); hold on
        plot(x, C2 * TOA_END ./ x, 'LineWidth', 2.2); hold on
        ylabel('count', 'FontSize', 12, 'FontName', 'cambria');
        xlabel('PRI / us', 'FontSize', 12, 'FontName', 'cambria');
        legend('actual count', 'first threshold', 'second threshold', ...
               'FontSize', 12, 'FontName', 'cambria');
        % judge whether exist valid pri
        nPRI = length(priRecord);
        isFindSeq = false;
        for k = 1 : nPRI
            onePRI = priRecord(k);
            if onePRI < minPRI, continue ; end
            if onePRI > maxPRI, break ; end
            if TOA_END / MIN_EVENTS < onePRI, break ; end
            [doublePRI, doubleK] = searchCloseTarget(2 * onePRI);
            if abs(doublePRI - 2 * onePRI) > TOLERANCE, continue ; end % there is no double pri
            oneEventCnt = indRecord(k).end - indRecord(k).begin + 1;
            doubleEventCnt = indRecord(doubleK).end - indRecord(doubleK).begin + 1;
            if oneEventCnt > C1 * TOA_END / onePRI && ...
               doubleEventCnt > C2 * TOA_END / doublePRI
                searchSeq(k);
                isFindSeq = true;
                C1 = C1 * DECAY;
                C2 = C2 * DECAY;
                break ;
            end
            
        end
        if isFindSeq
            diff = 1;
            record = [];
        else
            diff = diff + 1;
        end
    end

    function searchSeq(k)
       signal(end + 1).pri = priRecord(k);
       indArr = [];
       for j = indRecord(k).begin : indRecord(k).end
           indArr(end + 1) = record(j, 2); % first event index
           indArr(end + 1) = record(j, 3); % second event index
       end
       indArr = unique(indArr);
       indArr = sort(indArr);
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

    function [localPRI, localK] = searchCloseTarget(inArg)
        L = 1; R = length(priRecord);
        while R - L > 1
            M = ceil((L + R) / 2);
            if priRecord(M) < inArg
                L = M;
            elseif priRecord(M) > inArg
                R = M;
            else
                localPRI = priRecord(M);
                localK = M;
                return ;
            end
        end
        if inArg - priRecord(L) < priRecord(R) - inArg
            localPRI = priRecord(L);
            localK = L;
        else
            localPRI = priRecord(R);
            localK = R;
        end
    end
end
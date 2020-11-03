function priArr = sortingByTTP(toa, priRange, binNum)
% TTP algorithm
    c = 0.9;
    delta = 10;
    minPRI = priRange(1); maxPRI = priRange(2);
    nTOA = length(toa);
    towardMat = nan(nTOA, nTOA); inverseMat = nan(nTOA, nTOA);
    
    % toward planar transform 
    for i = 1 : nTOA
        ind = 1;
        for j = i + 1 : nTOA
            pri = toa(j) - toa(i);
            if pri >= minPRI && pri <= maxPRI
                towardMat(i, ind) = pri;
                ind = ind + 1;
            end
        end
    end
    
    % inverse planar transform
    for i = 1 : nTOA
        ind = 1;
        for j = i - 1 : -1 : 1
            pri = toa(i) - toa(j);
            if pri >= minPRI && pri <= maxPRI
                inverseMat(i, ind) = pri;
                ind = ind + 1;
            end
        end
    end
    
    towardMat = eraseNoiseByProb(towardMat);
    inverseMat = eraseNoiseByProb(inverseMat);
    
    % relative analysis
    for i = 1 : nTOA
        for j = 1 :nTOA
            u = towardMat(i, j);
            v = inverseMat(i, j);
            if ~isnan(u) && ~isnan(v)
                if abs(v - u) > delta
                    towardMat(i, j) = 0;
                end
            else
                break ;
            end
        end
    end
    
    % delete mirror pri
    priArr = zeros(1, nTOA);
    for i = 1 : nTOA
        row = towardMat(i, :);
        vals = row(~isnan(row));
        if ~isempty(vals)
            priArr(i) = vals(1);
        end
    end
    
    function mat = eraseNoiseByProb(mat)
        vals = mat(~isnan(mat));
        h = histogram(vals, binNum);
        meanVal = length(vals) / (maxPRI - minPRI);
        binWidth = h.BinWidth;
        for k = 1 : binNum
            binV = h.Values(k);
            if binV / binWidth <= c * meanVal
                mat(h.BinEdges(k) <= mat & h.BinEdges(k + 1) >= mat) = nan;
            end
        end
        close all;
    end
    
    
end
    



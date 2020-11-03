function [binVals, D] = sortingByPRITransform(toa, binNum, pri_range)
binVals = linspace(pri_range(1), pri_range(2), binNum + 1);
D = zeros(1, binNum);
nTOA = length(toa);

% enum every pulse interval and build histogram
for i = 1 : binNum
    for j = 1 : nTOA
        for k = j + 1 : nTOA
            pri = toa(k) - toa(j);
            if binVals(i) <= pri && binVals(i + 1) >= pri
                D(i) = D(i) + exp(2 * pi * 1j * toa(k) / (toa(k) - toa(j)));
            elseif pri > binVals(i + 1)
                break ;
            end
        end
    end
end
% bins middle value
binVals = binVals + (pri_range(2) - pri_range(1)) / binNum / 2;
binVals(end) = [];
D = abs(D);
end


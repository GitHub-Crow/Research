function arr = mergeSortedArray(arr1, arr2)
% merge arr1 and arr2 and two input array must be sorted
% example: arr1 = [1, 2, 10], arr2 = [3, 6, 14], then arr = [1, 2, 3, 6, 10, 14]
    n1 = length(arr1); n2 = length(arr2);
    k1 = 1; k2 = 1; k = 1;
    arr = nan(1, n1 + n2);
    
    while k1 <= n1 || k2 <= n2
        if k1 <= n1 && k2 <= n2
            if arr1(k1) < arr2(k2)
                arr(k) = arr1(k1);
                k1 = k1 + 1;
            else
                arr(k) = arr2(k2);
                k2 = k2 + 1;
            end
        elseif k1 <= n1
            arr(k) = arr1(k1);
            k1 = k1 + 1;
        else
            arr(k) = arr2(k2);
            k2 = k2 + 1;
        end
        k = k + 1;
    end
end


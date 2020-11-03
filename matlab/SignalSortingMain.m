clc; close all; clear ;
%% set parameters
TOA_RANGE = [0, 1e5];
PRI_DETECTED_RANGE  = [300, 700];

STABLE_PRI_1 = 400;
STABLE_PRI_2 = 510;

WOBBLE_PRI_1 = 400;
WOBBLE_PRI_2 = 610;
WOBBLE_RATE = 0.1;

STAGGER_PRI_1 = [400, 430, 450];
STAGGER_PRI_2 = [620, 640, 650];

SLIP_PRI_RANGE_1 = [400, 450];
SLIP_PRI_DELTA_1 = 5;
SLIP_PRI_RANGE_2 = [610, 670];
SLIP_PRI_DELTA_2 = 10;

SIN_PRI_PARAMS_1 = [400, 0.05, 50, 0]; % middle val, amplitude, cycle, init phase
SIN_PRI_PARAMS_2 = [610, 0.07, 70, 0];

%% generator TOA
stableTOA1 = generateTOA([STABLE_PRI_1], TOA_RANGE, PRIType.Stable);
stableTOA2 = generateTOA([STABLE_PRI_2], TOA_RANGE, PRIType.Stable);
stableTOA = mergeSortedArray(stableTOA1, stableTOA2);

wobbleTOA1 = generateTOA([WOBBLE_PRI_1, WOBBLE_RATE], TOA_RANGE, PRIType.Wobble);
wobbleTOA2 = generateTOA([WOBBLE_PRI_2, WOBBLE_RATE], TOA_RANGE, PRIType.Wobble);
wobbleTOA = mergeSortedArray(wobbleTOA1, wobbleTOA2);

staggerTOA1 = generateTOA([STAGGER_PRI_1], TOA_RANGE, PRIType.Stagger);
staggerTOA2 = generateTOA([STAGGER_PRI_2], TOA_RANGE, PRIType.Stagger);
staggerTOA = mergeSortedArray(staggerTOA1, staggerTOA2);

slipTOA1 = generateTOA([SLIP_PRI_RANGE_1, SLIP_PRI_DELTA_1], TOA_RANGE, PRIType.Slip);
slipTOA2 = generateTOA([SLIP_PRI_RANGE_2, SLIP_PRI_DELTA_2], TOA_RANGE, PRIType.Slip);
slipTOA = mergeSortedArray(slipTOA1, slipTOA2);

sinTOA1 = generateTOA([SIN_PRI_PARAMS_1], TOA_RANGE, PRIType.Sin);
sinTOA2 = generateTOA([SIN_PRI_PARAMS_2], TOA_RANGE, PRIType.Sin);
sinTOA = mergeSortedArray(sinTOA1, sinTOA2);

%% estimate PRI
ALGO_TYPE = 2;
% 0 for TTP transform; 1 for PRI transform; 2 for CDIF; 

if ALGO_TYPE == 0
    binNum = 10;
    stableD = sortingByTTP(stableTOA, PRI_DETECTED_RANGE, binNum);
    wobbleD = sortingByTTP(wobbleTOA, PRI_DETECTED_RANGE, binNum);
    staggerD = sortingByTTP(staggerTOA, PRI_DETECTED_RANGE, binNum);
    slipD = sortingByTTP(slipTOA, PRI_DETECTED_RANGE, binNum);
    sinD = sortingByTTP(sinTOA, PRI_DETECTED_RANGE, binNum);
    close all;
    figure ;
    cax = subplot(231);
    plotTOA_PRI(cax, stableTOA, stableD, 'stable PRI');
    ylim(PRI_DETECTED_RANGE);
    cax = subplot(232);
    plotTOA_PRI(cax, wobbleTOA, wobbleD, 'wobble PRI');
    ylim(PRI_DETECTED_RANGE);
    cax = subplot(233);
    plotTOA_PRI(cax, staggerTOA, staggerD, 'stagger PRI');
    ylim(PRI_DETECTED_RANGE);
    cax = subplot(234);
    plotTOA_PRI(cax, slipTOA, slipD, 'slip PRI');
    ylim(PRI_DETECTED_RANGE);
    cax = subplot(235);
    plotTOA_PRI(cax, sinTOA, sinD, 'sine PRI');
    ylim(PRI_DETECTED_RANGE);
    
elseif ALGO_TYPE == 1
    binNum = 200;
    [stableBins, stableD] = sortingByPRITransform(stableTOA, binNum, PRI_DETECTED_RANGE);
    [wobbleBins, wobbleD] = sortingByPRITransform(wobbleTOA, binNum, PRI_DETECTED_RANGE);
    [staggerBins, staggerD] = sortingByPRITransform(staggerTOA, binNum, PRI_DETECTED_RANGE);
    [slipBins, slipD] = sortingByPRITransform(slipTOA, binNum, PRI_DETECTED_RANGE);
    [sinBins, sinD] = sortingByPRITransform(sinTOA, binNum, PRI_DETECTED_RANGE);
    close all;
    figure ;
    cax = subplot(231);
    plotPRI_D(cax, stableBins, stableD, 'stable PRI'); 
    cax = subplot(232);
    plotPRI_D(cax, wobbleBins, wobbleD, 'wobble PRI');
    cax = subplot(233);
    plotPRI_D(cax, staggerBins, staggerD, 'stagger PRI');
    cax = subplot(234);
    plotPRI_D(cax, slipBins, slipD, 'slip PRI');
    cax = subplot(235);
    plotPRI_D(cax, sinBins, sinD, 'sine PRI');
    
elseif ALGO_TYPE == 2
    stableSig = sortingByCDIF(stableTOA, PRI_DETECTED_RANGE);
    wobbleSig = sortingByCDIF(wobbleTOA, PRI_DETECTED_RANGE);
    staggerSig = sortingByCDIF(staggerTOA, PRI_DETECTED_RANGE);
    slipSig = sortingByCDIF(slipTOA, PRI_DETECTED_RANGE);
    sinSig = sortingByCDIF(sinTOA, PRI_DETECTED_RANGE);
    
    close all;
    figure ;
    cax = subplot(231);
    plotTOA(cax, stableSig, 'stable PRI');
    cax = subplot(232);
    plotTOA(cax, wobbleSig, 'wobble PRI');
    cax = subplot(233);
    plotTOA(cax, staggerSig, 'stagger PRI');
    cax = subplot(234);
    plotTOA(cax, slipSig, 'slip PRI');
    cax = subplot(235);
    plotTOA(cax, sinSig, 'sine PRI');
    
end

%% plot function
function plotTOA(cax, sig, txt)
    axes(cax);
    nSig = length(sig);
    colorList = rand(nSig, 3);
    priList = strings;
    priList(end) = [];
    for k = 1 : nSig
        nSeq = length(sig(k).seq);
        stem(sig(k).seq, ones(1, nSeq), ...
             'MarkerFaceColor', colorList(k, :), ...
             'MarkerEdgeColor', 'k', ...
             'Marker', 'diamond'); hold on
        priList(end + 1) = num2str(sig(k).pri);
    end
    ylim([0, 2]);
    title(txt, 'FontSize', 14, 'FontName', 'cambria');
    ylabel('Value', 'FontSize', 12, 'FontName', 'cambria');
    xlabel('TOA / us', 'FontSize', 12, 'FontName', 'cambria');
    legend(priList, 'FontSize', 8, 'FontName', 'cambria');
end

function plotTOA_PRI(cax, x, y, txt)
    axes(cax);
    scatter(x, y, 'o', ...
            'MarkerEdgeColor', [0 .5 .5],...
            'MarkerFaceColor', [0 .7 .7],...
            'LineWidth', 1); 
    title(txt, 'FontSize', 14, 'FontName', 'cambria');
    ylabel('PRI / us', 'FontSize', 12, 'FontName', 'cambria');
    xlabel('TOA / us', 'FontSize', 12, 'FontName', 'cambria');
end

function plotPRI_D(cax, x, y, txt)
    axes(cax);
    plot(x, y, '-r',...
        'LineWidth', 2,...
        'MarkerSize', 1,...
        'MarkerEdgeColor', 'k',...
        'MarkerFaceColor', 'y');
    title(txt, 'FontSize', 14, 'FontName', 'cambria');
    ylabel('D', 'FontSize', 12, 'FontName', 'cambria');
    xlabel('PRI / us', 'FontSize', 12, 'FontName', 'cambria');
end

clc; clear ; close all;
%% set parameters
T = 4e-6;
FS = 100e6;
% single carrier frequency signal
PARAM_CW.amp = 1;
PARAM_CW.freq = 30e6;
PARAM_CW.phase = 0;
PARAM_CW.T = T;
PARAM_CW.FS = FS;
PARAM_CW.type = ModSigGen.CW;

% linear frequency modulation signal
PARAM_LFM.amp = 1;
PARAM_LFM.freq = 30e6;
PARAM_LFM.phase = 0;
PARAM_LFM.rate = 1.5e12;
PARAM_LFM.T = T;
PARAM_LFM.FS = FS;
PARAM_LFM.type = ModSigGen.LFM;

% phase coded signal
PARAM_PSK.amp = 1;
PARAM_PSK.freq = 30e6;
PARAM_PSK.phase = 0;
PARAM_PSK.code = [1, 1, 1, 0, 0, 1, 0, 1];
PARAM_PSK.encodingNum = 2;
PARAM_PSK.T = T;
PARAM_PSK.FS = FS;
PARAM_PSK.type = ModSigGen.PSK;

% freqency coded signal
PARAM_FSK.amp = 1;
PARAM_FSK.freqList = [20e6, 30e6, 25e6, 35e6];
PARAM_FSK.phase = 0;
PARAM_FSK.T = T;
PARAM_FSK.FS = FS;
PARAM_FSK.type = ModSigGen.FSK;

% frequency diversion signal
PARAM_FD.amp = 1;
PARAM_FD.freqList = [25e6, 30e6, 35e6];
PARAM_FD.phase = 0;
PARAM_FD.T = T;
PARAM_FD.FS = FS;
PARAM_FD.type = ModSigGen.FD;

%% generate modulation signal
PARAM_CW.S = ModSigGen.generate(PARAM_CW);
PARAM_LFM.S = ModSigGen.generate(PARAM_LFM);
PARAM_PSK.S = ModSigGen.generate(PARAM_PSK);
PARAM_FSK.S = ModSigGen.generate(PARAM_FSK);
PARAM_FD.S = ModSigGen.generate(PARAM_FD);

%% plot 
figure ;
cax = subplot(521);
plotTDS(PARAM_CW, cax);
cax = subplot(522);
plotFDS(PARAM_CW, cax);

cax = subplot(523);
plotTDS(PARAM_LFM, cax);
cax = subplot(524);
plotFDS(PARAM_LFM, cax);

cax = subplot(525);
plotTDS(PARAM_FSK, cax);
cax = subplot(526);
plotFDS(PARAM_FSK, cax);

cax = subplot(527);
plotTDS(PARAM_FD, cax);
cax = subplot(528);
plotFDS(PARAM_FD, cax);

cax = subplot(529);
plotTDS(PARAM_PSK, cax);
cax = subplot(5,2,10);
plotFDS(PARAM_PSK, cax);


function plotFDS(param, ax)
    f = linspace(-param.FS / 2, param.FS / 2, param.T * param.FS);
    axes(ax);
    y = abs(fftshift(fft(param.S)));
    plot(1e-6 .* f, y, '-r', 'LineWidth', 2);
    title(strcat(getTypeName(param.type), ' FD Signal'),...
         'FontSize', 14, 'FontName', 'cambria');
    ylabel('Amp / V', 'FontSize', 12, 'FontName', 'cambria');
    xlabel('f / MHz', 'FontSize', 12, 'FontName', 'cambria');
end

function plotTDS(param, ax)
    t = linspace(0, param.T, param.T * param.FS);
    axes(ax);
    plot(1e6 .* t, param.S, '-b', 'LineWidth', 2);
    title(strcat(getTypeName(param.type), ' TD Signal'),...
         'FontSize', 14, 'FontName', 'cambria');
    ylabel('PA / V', 'FontSize', 12, 'FontName', 'cambria');
    xlabel('t / us', 'FontSize', 12, 'FontName', 'cambria');
end
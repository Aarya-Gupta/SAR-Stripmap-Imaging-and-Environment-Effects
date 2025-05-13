

clear all; close all; clc;
tStart = tic;

disp('Initializing SAR simulation...');

%% 1. Define Radar System Parameters
speedOfLight = physconst('LightSpeed'); % Propagation speed (m/s)
carrierFreq = 5e9;                     % Center frequency: 5 GHz
wavelength = speedOfLight / carrierFreq; % Wavelength (m)
rangeRes = 2.5;                        % Range resolution: 2.5 m
azimuthRes = 2.5;                      % Azimuth resolution: 2.5 m
bandwidth = speedOfLight / (2 * rangeRes); % Signal bandwidth (Hz)
pulseFreq = 1200;                      % Pulse repetition frequency: 1200 Hz
pulseWidth = 2.5e-6;                   % Pulse duration: 2.5 µs
sampleRate = 150e6;                    % Sampling rate: 150 MHz
chirpRate = bandwidth / pulseWidth;    % Chirp rate (Hz/s)

%% 2. Configure Platform Motion
platformSpeed = 120;                   % Platform velocity: 120 m/s (y-axis)
altitude = 600;                        % Platform altitude: 600 m
flightTime = 1.8;                      % Flight duration: 1.8 s
radar = phased.Platform('InitialPosition', [0; -180; altitude], ...
                       'Velocity', [0; platformSpeed; 0]);

pulseInterval = 1 / pulseFreq;         % Time between pulses
numPulses = round(flightTime / pulseInterval); % Total pulses
maxDistance = 1200;                    % Maximum observable range: 1200 m
numSamples = ceil((2 * maxDistance / speedOfLight) * sampleRate); % Range samples
timeSamples = (0:numSamples-1) / sampleRate; % Fast-time vector

%% 3. Create Digital Elevation Model with Custom Terrain
gridSize = 80;                         % DEM grid size: 80x80 m
gridRes = 2;                           % DEM resolution: 2 m
[X, Y] = meshgrid(linspace(-gridSize/2, gridSize/2, gridSize/gridRes));

% Define Gaussian hills
hillCenters = [0, 0; 25, 25; -25, -25];    % [x, y] coordinates in meters
hillHeights = [50, 40, 60];                % Heights in meters
hillSigmas = [15, 12, 18];                 % Widths in meters
Z = zeros(size(X));                        % Initialize elevation grid
for k = 1:3
    X0 = hillCenters(k,1);                 % Hill center x-coordinate
    Y0 = hillCenters(k,2);                 % Hill center y-coordinate
    A = hillHeights(k);                    % Hill height
    sigma = hillSigmas(k);                 % Hill width
    Z = Z + A * exp(-((X - X0).^2 + (Y - Y0).^2) / (2 * sigma^2));
end
Z(Z < 0) = 0;                          % Ensure positive elevations

% SAR imaging grid
xGrid = X(1,:);                        % Cross-range axis
yGrid = Y(:,1);                        % Along-track axis
elevGrid = Z;                          % Elevation map

%% 4. Convert DEM to Scattering Points
scatterPoints = [];
scatterRCS = [];
sampleStep = 5;                        % Sample every 5th point (10 m spacing)
for row = 1:sampleStep:size(Z,1)
    for col = 1:sampleStep:size(Z,2)
        xPos = X(row,col);
        yPos = Y(row,col);
        zPos = Z(row,col);
        if abs(xPos) <= 40 && abs(yPos) <= 40
            scatterPoints = [scatterPoints, [xPos; yPos; zPos]];
            rcsVal = 0.2 + 0.8 * (zPos - min(Z(:))) / (max(Z(:)) - min(Z(:)));
            scatterRCS = [scatterRCS, rcsVal];
        end
    end
end
numScatterers = size(scatterPoints, 2);

%% 5. Initialize Radar Components
chirp = phased.LinearFMWaveform('SampleRate', sampleRate, ...
                               'PulseWidth', pulseWidth, ...
                               'PRF', pulseFreq, ...
                               'SweepBandwidth', bandwidth);
propChannel = phased.FreeSpace('PropagationSpeed', speedOfLight, ...
                              'OperatingFrequency', carrierFreq, ...
                              'SampleRate', sampleRate, ...
                              'TwoWayPropagation', true);
tx = phased.Transmitter('PeakPower', 1200, 'Gain', 45);
rx = phased.ReceiverPreamp('Gain', 45, 'NoiseFigure', 0);

%% 6. Simulate Radar Echoes
echoData = zeros(numSamples, numPulses);
disp(['Simulating ', num2str(numPulses), ' radar pulses...']);
parfor pulseIdx = 1:numPulses
    [platformPos, platformVel] = radar(pulseInterval);
    txPulse = chirp();
    txPulse = txPulse(1:numSamples);
    txSignal = tx(txPulse);
    echoTemp = zeros(numSamples, 1);
    for scatterIdx = 1:numScatterers
        [~, angles] = rangeangle(scatterPoints(:,scatterIdx), platformPos);
        propSignal = propChannel(txSignal, platformPos, scatterPoints(:,scatterIdx), ...
                                 platformVel, [0;0;0]);
        echoTemp = echoTemp + propSignal * sqrt(scatterRCS(scatterIdx));
    end
    echoData(:,pulseIdx) = rx(echoTemp);
end
disp(['Echo simulation completed in ', num2str(toc(tStart)), ' seconds']);

% Platform positions
platformY = -180 + (0:numPulses-1) * platformSpeed * pulseInterval;

%% 7. Apply Range Compression
disp('Performing range compression...');
rangeProcessor = phased.RangeResponse('RangeMethod', 'Matched filter', ...
                                     'PropagationSpeed', speedOfLight, ...
                                     'SampleRate', sampleRate);
filterCoeff = getMatchedFilter(chirp);
[compressedData, rangeGrid] = rangeProcessor(echoData, filterCoeff);

%% 8. Form SAR Image via Backprojection
disp('Generating SAR image...');
sarImage = backproject(compressedData, rangeGrid, platformY, xGrid, yGrid, ...
                       elevGrid, altitude, speedOfLight, carrierFreq);

%% 9. Display Results
% Plot DEM
figure;
surf(X, Y, Z, 'EdgeColor', 'none');
title('Digital Elevation Model');
xlabel('Cross-range (m)');
ylabel('Along-track (m)');
zlabel('Height (m)');
axis tight;
view(45, 50);
colorbar;

% Plot SAR Image
figure('Position', [200, 200, 1000, 600]);
imagesc(xGrid, yGrid, abs(sarImage));
axis xy tight;
title('SAR Stripmap Image');
xlabel('Cross-range (m)');
ylabel('Along-track (m)');
colorbar;

disp(['Simulation completed in ', num2str(toc(tStart)), ' seconds']);

%% Backprojection Function
function sarImage = backproject(data, ranges, platformPos, xCoords, yCoords, ...
                               elevMap, height, c, freq)
    [rows, cols] = size(elevMap);
    sarImage = zeros(rows, cols);
    for r = 1:rows
        for c = 1:cols
            x = xCoords(c);
            y = yCoords(r);
            z = elevMap(r,c);
            distances = sqrt(x^2 + (y - platformPos).^2 + (z - height)^2);
            delays = 2 * distances / c;
            for p = 1:length(platformPos)
                sample = interp1(ranges, data(:,p), distances(p), 'linear', 0);
                sarImage(r,c) = sarImage(r,c) + sample * exp(1j*2*pi*freq*delays(p));
            end
        end
    end
end
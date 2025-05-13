% Main script for SAR simulation and processing
function SAR_Simulation()
    % Step 1: Set up parameters
    params = setupParameters();
    
    % Step 2: Generate waveform
    waveform = generateWaveform(params);
    
    % Step 3: Initialize platform and target
    [radarPlatform, pointTargets] = initializePlatforms(params);
    
    % Step 4: Collect SAR data
    rxsig = collectSARData(params, waveform, radarPlatform, pointTargets);
    
    % Step 5: Perform pulse compression
    [cdata, rnggrid] = performPulseCompression(params, waveform, rxsig);
    
    % Step 6: Apply back-projection algorithm
    bpa_processed = applyBackProjection(params, cdata, rnggrid);
    
    % Step 7: Visualize the focused image
    visualizeImage(bpa_processed);
end

% Function to set up parameters
function params = setupParameters()
    params.c = physconst('LightSpeed');         % Speed of light
    params.fc = 6e9;                            % Center frequency (6 GHz)
    params.lambda = params.c / params.fc;       % Wavelength
    params.rangeResolution = 3;                 % Range resolution (m)
    params.crossRangeResolution = 3;            % Cross-range resolution (m)
    params.bw = params.c / (2 * params.rangeResolution); % Bandwidth
    params.prf = 3000;                          % Pulse Repetition Frequency (Hz)
    params.tpd = 3e-6;                          % Pulse duration (s)
    params.fs = 120e6;                          % Sampling rate (Hz)
    params.K = params.bw / params.tpd;          % Chirp rate
    params.v = 300;                             % Aircraft velocity (m/s)
    params.h = 1500;                            % Altitude (m)
    params.flightDuration = 4;                  % Flight duration (s)
    params.maxRange = 2500;                     % Maximum range (m)
    params.truncrangesamples = ceil((2 * params.maxRange / params.c) * params.fs);
    params.fastTime = (0 : 1/params.fs : (params.truncrangesamples - 1) / params.fs)';
    params.Rc = 1000;                           % Reference range (m)
    params.P_tx = 1500;                         % Transmit power (W)
    params.G = 60;                              % Antenna gain (dBi)
end

% Function to generate waveform
function waveform = generateWaveform(params)
    waveform = phased.LinearFMWaveform( ...
        'SampleRate', params.fs, ...
        'PulseWidth', params.tpd, ...
        'PRF', params.prf, ...
        'SweepBandwidth', params.bw, ...
        'SweepDirection', 'Up');
end

% Function to initialize platforms
function [radarPlatform, pointTargets] = initializePlatforms(params)
    radarPlatform = phased.Platform( ...
        'InitialPosition', [0; -200; params.h], ...
        'Velocity', [0; params.v; 0]);
    targetpos = [1000; 0; 0];                   % Target at [1000, 0, 0]
    targetvel = [0; 0; 0];                      % Stationary target
    pointTargets = phased.Platform('InitialPosition', targetpos, 'Velocity', targetvel);
end

% Function to collect SAR data
function rxsig = collectSARData(params, waveform, radarPlatform, pointTargets)
    slowTime = 1 / params.prf;
    numpulses = round(params.flightDuration / slowTime);
    rxsig = zeros(params.truncrangesamples, numpulses);
    refangle = 0;                               % 0° squint angle
    
    % Define components
    channel = phased.FreeSpace( ...
        'PropagationSpeed', params.c, ...
        'OperatingFrequency', params.fc, ...
        'SampleRate', params.fs, ...
        'TwoWayPropagation', true);
    target = phased.RadarTarget('OperatingFrequency', params.fc, 'MeanRCS', 1);
    radiator = phased.Radiator('Sensor', phased.IsotropicAntennaElement, ...
        'OperatingFrequency', params.fc);
    collector = phased.Collector('Sensor', phased.IsotropicAntennaElement, ...
        'OperatingFrequency', params.fc);
    transmitter = phased.Transmitter('PeakPower', params.P_tx, 'Gain', params.G);
    receiver = phased.ReceiverPreamp('Gain', params.G, 'NoiseFigure', 0);
    
    % Data collection loop
    for ii = 1:numpulses
        [radarpos, radarvel] = radarPlatform(slowTime);
        [targetpos, targetvel] = pointTargets(slowTime);
        pulse = waveform(); pulse = pulse(1:params.truncrangesamples);
        txsig = transmitter(pulse);
        txsig = radiator(txsig, refangle);
        txsig = channel(txsig, radarpos, targetpos, radarvel, targetvel);
        txsig = target(txsig);
        txsig = collector(txsig, refangle);
        rxsig(:, ii) = receiver(txsig);
    end
end

% Function to perform pulse compression
function [cdata, rnggrid] = performPulseCompression(params, waveform, rxsig)
    pulseCompression = phased.RangeResponse( ...
        'RangeMethod', 'Matched filter', ...
        'PropagationSpeed', params.c, ...
        'SampleRate', params.fs);
    mfCoeff = getMatchedFilter(waveform);
    [cdata, rnggrid] = pulseCompression(rxsig, mfCoeff);
end

% Function to apply back-projection algorithm
function bpa_processed = applyBackProjection(params, sigdata, rnggrid)
    fastTime = params.fastTime;
    fc = params.fc;
    fs = params.fs;
    prf = params.prf;
    v = params.v;
    crossRangeResolution = params.crossRangeResolution;
    c = params.c;
    
    bpa_processed = helperBackProjection(sigdata, rnggrid, fastTime, fc, fs, prf, v, crossRangeResolution, c);
end

% Helper function for Back-Projection Algorithm
function data = helperBackProjection(sigdata, rnggrid, fastTime, fc, fs, prf, speed, crossRangeResolution, c)
    data = zeros(size(sigdata));
    azimuthDist = -200:speed/prf:200;           % Azimuth positions
    rangelims = [700 1400];                     % Range window (m)
    crossrangelims = [-10 10];                  % Cross-range window (m)
    
    rangeIdx = [find(rnggrid > rangelims(1), 1), find(rnggrid < rangelims(2), 1, 'last')];
    crossrangeIdxStart = find(azimuthDist > crossrangelims(1), 1);
    crossrangeIdxStop = find(azimuthDist < crossrangelims(2), 1, 'last');
    
    for i = rangeIdx(1):rangeIdx(2)
        posy = c * fastTime(i) / 2;             % Range from time
        lsynth = (c / fc) * posy / (2 * crossRangeResolution); % Synthetic aperture length
        lsar = round(lsynth * length(azimuthDist) / abs(azimuthDist(end)));
        lsar = lsar + mod(lsar, 2);            % Ensure even length
        hn = hanning(lsar).';                   % Hanning window
        
        for j = crossrangeIdxStart:crossrangeIdxStop
            posx = azimuthDist(j);
            count = 0;
            
            for k = j - lsar/2 + 1 : j + lsar/2
                if k < 1 || k > size(sigdata, 2)
                    continue;
                end
                td = 2 * sqrt((azimuthDist(k) - posx)^2 + posy^2) / c; % Two-way delay
                cell = round(td * fs) + 1;
                
                if cell >= 1 && cell <= size(sigdata, 1)
                    signal = sigdata(cell, k);
                    count = count + hn(k - (j - lsar/2)) * signal * exp(1j * 2 * pi * fc * td);
                end
            end
            data(i, j) = count;
        end
    end
end

% Function to visualize the focused image
function visualizeImage(bpa_processed)
    figure;
    imagesc(abs(bpa_processed(600:1400, 1700:2300)));
    title('SAR Data focused using Back-Projection Algorithm');
    xlabel('Cross-Range Samples');
    ylabel('Range Samples');
end
% Original parameters (unchanged)
c = physconst('LightSpeed');
lambda = c / 6e9;         % Wavelength for 6 GHz radar
fc = 6e9;                 % Center frequency
rangeResolution = 3;      % in meters
crossRangeResolution = 3;
bw = c / (2 * rangeResolution);  % Bandwidth
prf = 3000;               % Pulse Repetition Frequency
tpd = 3e-6;               % Pulse duration
fs = 120e6;               % Sampling rate
K = bw / tpd;             % Chirp rate

% Waveform
waveform = phased.LinearFMWaveform( ...
    'SampleRate', fs, ...
    'PulseWidth', tpd, ...
    'PRF', prf, ...
    'SweepBandwidth', bw, ...
    'SweepDirection', 'Up');

% Platform parameters
v = 300;                   % Aircraft velocity in m/s
h = 1500;                  % Altitude in meters
flightDuration = 4;        % seconds
radarPlatform = phased.Platform( ...
    'InitialPosition', [0; -200; h], ...
    'Velocity', [0; v; 0]);

slowTime = 1 / prf;
numpulses = round(flightDuration / slowTime);
maxRange = 2500;
truncrangesamples = ceil((2 * maxRange / c) * fs);
fastTime = (0 : 1/fs : (truncrangesamples - 1) / fs);
Rc = 1000;                % Reference range

% Free-space channel
channel = phased.FreeSpace( ...
    'PropagationSpeed', c, ...
    'OperatingFrequency', fc, ...
    'SampleRate', fs, ...
    'TwoWayPropagation', true);

% Single Target Setup
targetpos = [1000; 0; 0];  % Single target at [1000, 0, 0]
targetvel = [0; 0; 0];     % Stationary target
target = phased.RadarTarget('OperatingFrequency', fc, 'MeanRCS', 1);
pointTargets = phased.Platform('InitialPosition', targetpos, 'Velocity', targetvel);

% Antenna gains and power
P_tx = 1500; % Transmit power in watts
G = 60;      % Antenna gain in dBi
radiator = phased.Radiator('Sensor', phased.IsotropicAntennaElement, ...
    'OperatingFrequency', fc);
collector = phased.Collector('Sensor', phased.IsotropicAntennaElement, ...
    'OperatingFrequency', fc);
transmitter = phased.Transmitter('PeakPower', P_tx, 'Gain', G);
receiver = phased.ReceiverPreamp('Gain', G, 'NoiseFigure', 0);

% Pulse Compression setup
pulseCompression = phased.RangeResponse( ...
    'RangeMethod', 'Matched filter', ...
    'PropagationSpeed', c, ...
    'SampleRate', fs);
mfCoeff = getMatchedFilter(waveform);

% Environmental conditions to simulate
conditions = {'none', 'rain', 'fog', 'snowfall'};
rain_rate = 10;          % mm/hr for rain
visibility = 1000;       % meters for fog
snowfall_rate = 5;       % mm/hr equivalent for snowfall

% Store processed data for each condition
bpa_processed_all = cell(length(conditions), 1);

% Loop over environmental conditions
for cond_idx = 1:length(conditions)
    env_condition = conditions{cond_idx};
    
    % SAR data collection with environmental effects
    rxsig = zeros(truncrangesamples, numpulses);
    refangle = 0;  % 0° squint for single target
    
    for ii = 1:numpulses
        [radarpos, radarvel] = radarPlatform(slowTime);
        [targetpos, targetvel] = pointTargets(slowTime);
        
        [targetRange, targetAngle] = rangeangle(targetpos, radarpos);
        
        pulse = waveform(); pulse = pulse(1:truncrangesamples);
        txsig = transmitter(pulse);
        txsig = radiator(txsig, refangle);
        
        % Propagate signal through channel
        txsig = channel(txsig, radarpos, targetpos, radarvel, targetvel);
        
        % Apply environmental attenuation
        if strcmp(env_condition, 'rain')
            attenuation = calculate_rain_attenuation(targetRange, rain_rate, fc);
        elseif strcmp(env_condition, 'fog')
            attenuation = calculate_fog_attenuation(targetRange, visibility, fc);
        elseif strcmp(env_condition, 'snowfall')
            attenuation = calculate_snowfall_attenuation(targetRange, snowfall_rate, fc);
        else
            attenuation = 1;  % No attenuation for 'none'
        end
        
        % Apply two-way attenuation to the signal
        txsig = txsig * attenuation;
        txsig = target(txsig);
        txsig = collector(txsig, refangle);
        rxsig(:, ii) = receiver(txsig);
    end
    
    % Pulse Compression
    [cdata, rnggrid] = pulseCompression(rxsig, mfCoeff);
    
    % Back-Projection processing
    bpa_processed = helperBackProjection(cdata, rnggrid, fastTime, fc, fs, prf, v, crossRangeResolution, c);
    bpa_processed_all{cond_idx} = bpa_processed;
    
    % Plot the SAR image for this condition
    figure;
    imagesc(abs(bpa_processed(600:1400, 1700:2300)));
    title(['SAR Image - ', env_condition, ' Condition']);
    xlabel('Cross-Range Samples'); ylabel('Range Samples');
    colormap('jet'); colorbar;
end

% Helper functions for environmental attenuation
function atten = calculate_rain_attenuation(range, rain_rate, fc)
    % Simple rain attenuation model (e.g., based on ITU-R P.838)
    alpha = 0.01 * rain_rate^0.6;  % Specific attenuation (dB/km)
    atten_db = 2 * alpha * (range / 1000);  % Two-way attenuation in dB
    atten = 10^(-atten_db / 10);  % Convert to linear scale
end

function atten = calculate_fog_attenuation(range, visibility, fc)
    % Simple fog attenuation model based on visibility
    if visibility > 1000
        atten = 1;  % Negligible attenuation
    else
        alpha = 0.1 / visibility;  % Rough estimate (dB/km)
        atten_db = 2 * alpha * (range / 1000);  % Two-way attenuation in dB
        atten = 10^(-atten_db / 10);
    end
end

function atten = calculate_snowfall_attenuation(range, snowfall_rate, fc)
    % Simple snowfall attenuation model
    alpha = 0.005 * snowfall_rate^0.7;  % Specific attenuation (dB/km)
    atten_db = 2 * alpha * (range / 1000);  % Two-way attenuation in dB
    atten = 10^(-atten_db / 10);  % Convert to linear scale
end

% Helper function for Back-Projection Algorithm (unchanged)
function data = helperBackProjection(sigdata, rnggrid, fastTime, fc, fs, prf, speed, crossRangeResolution, c)
    data = zeros(size(sigdata));
    azimuthDist = -200:speed/prf:200;

    rangelims = [700 1400];
    crossrangelims = [-10 10];

    rangeIdx = [find(rnggrid > rangelims(1), 1), find(rnggrid < rangelims(2), 1, 'last')];
    crossrangeIdxStart = find(azimuthDist > crossrangelims(1), 1);
    crossrangeIdxStop = find(azimuthDist < crossrangelims(2), 1, 'last');

    for i = rangeIdx(1):rangeIdx(2)
        posy = c * fastTime(i) / 2;
        lsynth = (c / fc) * posy / (2 * crossRangeResolution);
        lsar = round(lsynth * length(azimuthDist) / abs(azimuthDist(end)));
        lsar = lsar + mod(lsar, 2);
        hn = hanning(lsar).';

        for j = crossrangeIdxStart:crossrangeIdxStop
            posx = azimuthDist(j);
            count = 0;

            for k = j - lsar/2 + 1 : j + lsar/2
                if k < 1 || k > size(sigdata, 2)
                    continue;
                end

                td = 2 * sqrt((azimuthDist(k) - posx)^2 + posy^2) / c;
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
% Define parameters
fs = 1000;          % Sampling frequency (Hz)
t = 0:1/fs:1;       % Time vector (0 to 1 second)
freq1 = 5;          % Frequency of the first ramp waveform (Hz)
freq2 = 5;          % Frequency of the second ramp waveform (Hz)
delay = 0.2;        % Time delay between the two waveforms (seconds)

% Generate the first ramp waveform
waveform1 = sawtooth(2*pi*freq1*t);

% Generate the second ramp waveform with time delay
waveform2 = sawtooth(2*pi*freq2*(t - delay));

% Calculate the instantaneous frequency of both waveforms
instantaneous_freq1 = diff(unwrap(angle(hilbert(waveform1)))) / (2*pi);
instantaneous_freq2 = diff(unwrap(angle(hilbert(waveform2)))) / (2*pi);

% Calculate the frequency difference
frequency_difference = instantaneous_freq2 - instantaneous_freq1;

% Plot the results
subplot(3,1,1);
plot(t, waveform1);
title('Waveform 1');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3,1,2);
plot(t, waveform2);
title('Waveform 2 (Time-Delayed)');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3,1,3);
plot(t(1:end-1), frequency_difference);
title('Frequency Difference');
xlabel('Time (s)');
ylabel('Frequency Difference (Hz)');

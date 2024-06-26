% Parameters %% HowMLGeron can be employed
uav_speed = 100; % UAV speed in m/s
angle_of_arrival = 30; % Angle of arrival in degrees

fc = 28e9;   %2.4e9; % Carrier frequencyHz%[0.6 2.4 4.2 28]*1e9

c = 3e8;    % Speed of light in m/s

lambda = c / fc; % Wavelength

num_elements_x = 4; % Number of elements along x-axis

num_elements_y = 4; % Number of elements along y-axis

element_spacing_x = lambda / 2; % Half-wavelength spacing along x-axis

element_spacing_y = lambda / 2; % Half-wavelength spacing along y-axis

fs = 8e9; % Sampling frequency

duration = 1; % Signal duration in seconds

SNR = 20; % Signal-to-noise ratio in dB

% UAV speed parameters
true_speed = 100; % True speed of the UAV in m/s

doppler_shift = uav_speed / c;

baseband_frequency = fc - doppler_shift;
%%
chirp_bandwidth = 1e6; % Chirp bandwidth in Hz
% t = 0:1/1e3:2;
% y = chirp(t,100,1,200,'quadratic');
baseband_GHz = baseband_frequency/1e9;
bandwidth_KHz = chirp_bandwidth/1e3;
chirp_signal_baseband = ...
chirp(0:1e8/fs:duration, baseband_GHz, duration, bandwidth_KHz,'quadratic');

% Create a phased array
array_position_x = (0:(num_elements_x - 1)) * element_spacing_x;
array_position_y = (0:(num_elements_y - 1)) * element_spacing_y;

% Calculate Doppler shift for each element
doppler_shift_elements = zeros(num_elements_x, num_elements_y);
for m = 2:num_elements_x
    for n = 2:num_elements_y
        angle_x = sind(angle_of_arrival) * (m - 1) * element_spacing_x / sqrt((m - 1) * element_spacing_x)^2;
        angle_y = sind(angle_of_arrival) * (n - 1) * element_spacing_y / sqrt((n - 1) * element_spacing_y)^2;
        
        relative_velocity = uav_speed * (cosd(angle_x) + cosd(angle_y));
        doppler_shift_elements(m, n) = relative_velocity / c;
    end
end

% Generate received signals at each element
received_signals = zeros(num_elements_x, num_elements_y, length(chirp_signal_baseband));
for m = 2:num_elements_x
    for n = 2:num_elements_y
        received_signals(m, n, :) = chirp_signal_baseband .* exp(1i * 2 * pi * doppler_shift_elements(m, n) * (0:1e8/fs:duration));
    end
end

% Display results
disp(['UAV Speed: ' num2str(uav_speed) ' m/s']);
disp('Doppler Shift for Each Element:');
disp(doppler_shift_elements);

% Plot results
figure;
subplot(2, 1, 1);
plot(0:1e8/fs:duration, real(chirp_signal_baseband));
title('Transmitted Chirp Signal');
xlabel('Time (s)');
ylabel('Amplitude');

dim = num_elements_x;
subplot(2, 1, 2);
stem3(1:num_elements_x, 1:num_elements_y, doppler_shift_elements(:,dim)/1e-9);
title('Doppler Shift for Each Array Element');
xlabel('Element Index (x)');
ylabel('Element Index (y)');
zlabel('Doppler Shift (Hz)');


% In this code, a linear chirp signal is generated and transmitted from a moving UAV. 
% The received signals at each element of the phased array are calculated, and noise 
% is added to simulate a realistic scenario. Doppler estimation is performed using 
% the phased array and Doppler estimator, and the estimated UAV speed is then 
% compared with the true speed.





%% Test Phase Pulse compression
% For 1:NSAs, computeO/p of SAs for steerEachSA with phase comp for targ from each element of SA
%%%%%%

theta_degree= 45; %30;   %angle_of_arrival
v_UAV = 100;        % UAV velocity in m/s
  f_Hz = 28e9; %2.4e9; 1e9
for ai= 1: NSAs    

    array_x = arraystack_x (atemp, 1:Nx); %use atemp
    array_y = arraystack_y (atemp, 1:Ny); 
    steering_phase = SA_steer( array_x, array_y, -u_steer, -v_steer, k );
    
    atemp  = 10*ai + [1:1:Nx]; %next loop value
    tempSA = zeros (length(t), 1);

    for i= 1 : Nx  % Nx=Ny=10 %1st SA
    for j= 1 : Ny
    %[xQ yQ zQ] = fAzel2xy(azQ, elQ, r); %1.object position
    xQ = targX;    yQ = targY;    zQ = targZ;
   
    dx =xQ - array_x (i,j) ;%2. distance bw object and every element
    dy =yQ - array_y (i,j);
    dz =zQ;
    d = sqrt(dx.^2 + dy.^2 + dz.^2 ); %

    phase(i,j) = mod( (2*pi/lambda)* d, 2 * pi );%3.1. phase_shift based on pathshift
    phi_init_rad = phase(i,j);%3.2 chirp for this phase ip
    (i,j)
    %complex valued analog signal using the I/Q sg at the pulse compressor input
   %steering_vector.Direction = [angle_of_arrival; 0];
  env_norm = @(lfm_scal, phi_init_rad, t, mu, wc, SNR)...
        (SNR * (exp ( 1j *(wc * t + mu * (t .^ 2 ) / 2 + phi_init_rad )) ) / lfm_scal ); 
lfm_comp_env_norm = env_norm (lfm_scal,...
      phi_init_rad, t, mu, wc, SNR );
   
  %%% DOppler Shift
  
  doppler_shift = v_UAV*sin(theta)*f_Hz/ c;
  doppler_beta= exp(1i * 2 * pi * doppler_shift * (i+j));
 Doppler_lfm_comp_env_norm(i,j) = ...
 doppler_beta*lfm_comp_env_norm ;
    
   %end
    %%%%%
    I_lfm = real (lfm_comp_env_norm);
    Q_lfm = imag (lfm_comp_env_norm);
    chirp = I_lfm + 1i * Q_lfm;
    chirp = reshape( chirp, Nsam, 1 );
    noise = sqrt(2) * (randn(L,1) + 1j*randn(L,1) );
    chirp = chirp + noise; 
    
    Ep = 10 * log10( sum( lfm_comp_env_norm .* conj( lfm_comp_env_norm ) ) );%Energy
    x_signalstack(ai,i,j,:) = chirp  ; %NSAs x Nx x Ny x L

%steering_phase = SA_steer( array_x, array_y, -u_steer, -v_steer, k);SA_steer = @(x,y,u,v,k)(exp(-1j*k*(x.*u + y.*v)));
    tempSA = tempSA + steering_phase(i,j) * chirp;   %  coherently combine the elements of each subarray (beamforming)
    end
end


%%%%%%%%%%%%%
% Plot the array pattern
figure;
pattern(array_response, fc, -180:180, 0, 'PropagationSpeed', c, 'Type', 'powerdb', 'ElementSpacing', array_spacing);
title('Phased Array Antenna Pattern');

% Display the received signal
figure;
stem(1:num_elements, abs(received_signal));
xlabel('Antenna Element');
ylabel('Received Signal Magnitude');
title('Received Signal at Phased Array Elements');

% Visualize the array geometry
figure;
viewArray(array_response, 'ShowIndex', 'All', 'ShowNormals', true);
title('Phased Array Geometry');
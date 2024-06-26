close all
clear all
clc
%%
tic
fs          = 8e9;          %  sampling frequency
delF        = 1e9;
fc          = 0.001;         %  baseband frequency
wc          = 2 * pi * fc;
Ts          = 1 / fs; %  sampling period
T           = 1e-6;  %10e-6;           %  chirp duration
t           = 0 : Ts : T;
Nsam        = length( t );

SNR_dB      = 60;
phi_init_rad= 0;
BT          = delF * T;          %  time-bandwidth product
SNR         = 10 ^ ( SNR_dB / 20 ); %snr in voltage ratio
lfm_scal    = 1;









%  Calculate frequency sweep rate (radians/sec)
mu          = 2 * pi * delF / T;
%% Beam of beam of subarrays 1. the kth element of array belongs to the sth subarray the ith lement
% Nx0= 6 %30;
% Ny0= 6 %30;
% Nelem= Nx0 * Ny0;
% NSAs =9 %9;
%%
Nx0         = 2;
Ny0         = 2;
Nelem       = Nx0*Ny0;
NSAs        = 1
% Assuming 1st SA at center (0,0) (position in matrix form is 1st= 11)
Nx          = sqrt (Nelem /NSAs);       %10
Ny          = sqrt (Nelem/ NSAs);       %10
Nelem_SA    = Nx * Ny;
%%
% f_Hz      =1e9;
f_Hz        = 2.8e9;  
c           = 3e8; 
lambda      = c/f_Hz; 
%k= 2*pi/lambda; 
%spacing = lambda/2; % warning('off','all');

freq_sam    = 1351; %vary sample size and see impact and think
% Nelem     = 36; NelemOver2 = ( Nelem - 1 ) / 2; 
% NelemTot  = Nelem * Nelem;
BandWidth_Hz = 0;  %10e9;  %1e9;

fLO         = 2.65e9;
fHI         = fLO + BandWidth_Hz;
fMID        = 0.5 * ( fLO + fHI );
%Input range
f_Hz        = linspace( fLO, fHI, freq_sam ); 

lambdaHI    = c / fHI; 
spacing     = lambdaHI / 2;
lambdaMID   = c / fMID; 
k           = 2 * pi / lambdaHI;
save_array  =1;
%% clear tempx      %square NSAs_x = NSAs_y ; Nx = Ny
NSAs_x      = Nx0/Nx;     NSAs_y = NSAs/ NSAs_x;
a_x         = zeros(NSAs_x* Nx, Nx);
a_y         = zeros(NSAs_y* Ny, Ny);

temp        = [1:Nx];           %initial index, it starts from 1 not 0
g           = [0:Nx-1];            % initial elemt position in a SA, starts from 0

%create SA of Nx size and put values in orderexpand till the entire all SA
for ia_x    = 1: NSAs_x    % combine SA matx from 1 to NSAs_x
% %a_x( temp, 1:Nx ) = kron((grid),ones (Nx,1)) ;      %only a_x

[a_x(temp, 1:Nx) a_y(temp, 1:Ny)] = meshgrid(g , g) ;

temp        = Nx + temp  ;    %next loop value %   temp_x = Nx * ia_x + [1:Nx];
    g       = Nx + g;
end

array_x     = a_x * spacing; %30x10 %1st 2nd 3rd x/y of SA, rwise
array_y     = a_y * spacing ;

%% stack o/ps % stack periodic matrix row wise   % perm = repmat(eye(NSAs_x), [NSAs_x 1]);
%Arrange ones in SAs and expand till array
perm           = repmat( eye(NSAs_x * Nx), [NSAs_x 1] );
stack_a_x      =   perm * a_x;         %90x10 %all SA, rwise
arraystack_x   = stack_a_x * spacing;

stack_a_y      = perm * a_y;
arraystack_y   = stack_a_y * spacing;

%% plot all the SA in x-y plane
 plot_array    = 1;% Control statemnt
%plot_array    = 0;
if plot_array
for short = 1:1
% To plot entire array    
%origin_vec =  [0: Nx-1];
% origin_pos =  [-1 0 0; 0 0 0; 0 0 0]; %ehen Nx0 is multiple of 3
% origin_next_pos = eye (NSAs_x) +  origin_pos;
disp('Discreate blueprint of array and other arrays')
origin_vec =  [0: Nx-1];
origin_next_pos = ones(NSAs);
origin_pos     = 0;
origin_next_pos(1,1)= origin_pos;
initial_index = Nx * origin_next_pos * [0: NSAs_x-1 ].';
next_increment= repmat( initial_index, [1, Nx]  );

order_xy = repmat(origin_vec, [NSAs_x, 1])+ next_increment;
order_SA = order_xy* spacing;

[xall_SA yall_SA] = meshgrid( order_SA, order_SA);

%% plot max min cent0 of SA

atemp = [1:Nx];         %initial index, it starts from 1 not 0
ia    = NSAs ;
SA_x0 = zeros (NSAs,1) ;  
SA_y0 = zeros(NSAs,1);  array_x0= zeros(NSAs) ; array_y0 =zeros (NSAs);
%% To plot SA centers over the existing array
 for ai= 1: NSAs         % vary the SA from 1 to NSAs
    array_x0 = arraystack_x (atemp, 1:Nx); %use atemp 10x 10
    array_y0 = arraystack_y (atemp, 1:Ny); 
 
    SA_center= @(array_x)( min(array_x(:))+ ( max(array_x(:))-min(array_x(:)) )/2 );  
    SA_x0 (ai) = SA_center (array_x0) ;
    SA_y0 (ia) = SA_center (array_y0) ;     %y:oppositeOrder than x
        SA_x0         
     ia = NSAs - ai;
    atemp = Nx + atemp;         %10*ai + [1:1:Nx]; %next loop value
 end

%%
[SA_x00 SA_y00] = meshgrid (SA_x0, SA_y0);
figure(9)
plot( SA_x00, SA_y00, 'vm'); grid on
hold on
plot(xall_SA, yall_SA, 'ob'); grid on

% targR = 2e3;
% targX = 14e3;
% targY = 4e3;
end
end
if save_array
save uav_data_file_array.mat xall_SA yall_SA arraystack_x arraystack_y SA_x00 SA_y00 Nx0 Ny0 Nelem NSAs
else
display('Load the Array Layout')
load uav_data_file_array
end
%load uav_data_file_array
%% phase computation of target location& add possible Error for each element
L = length(t);
targY = 7e3/3.28084;          %[7e3 10.5e3 13e3 9.5e3 14.5e3]/3.28084;
targX = 3 * -14e3/3.28084;    %3 * [-14e3 -7e3 2e3 7e3 13e3]/3.28084;
targR = 200 * 1e2;            %[200 180 197.5 175 199.5] * 1e2;

%find TRUELocasn f trgt z-dim (m) from given range, x, and y (objectLocation)
targZ = sqrt(targR.^2 - targX.^2 - targY.^2);
targPhi = atan2(targY,targX);           % target angle phi from given x and y
targTh = acos(targZ./targR);            % target angle theta from given z and range

[ mpc_ux, mpc_vy ] = fThetaphi2uv( targTh, targPhi );
figure(29)
plot( targTh, targPhi, targZ, 'vm'); grid on
hold on
plot(mpc_vy, mpc_ux, 'ob'); grid on

%%
%UserDefFormula Error bw steeringDixn and trueSigDixn in 1dim(x) 
display('Steer the array')
Bfrac = 10;
%mainbeam_3dB = 0.886 * lambda / ( Nx0 * lambda );   
mainbeam_3dB =lambda / (Nelem*lambda);

% ux_aoa_err = ( mainbeam_3dB / Bfrac );      %rand( 1, array.nTarg );        %  user defined error between steering direction and true signal direction %( mainbeam_3dB / 8 )*[1 -1 1-1] in quadrant beam error
% vy_aoa_err = ( mainbeam_3dB / Bfrac );      %rand( 1, array.nTarg );
% initial_rms_error = sqrt( ux_aoa_err .^2 + vy_aoa_err .^2 );
% u_steer = mpc_ux + ux_aoa_err; %recv beam steering direction
% v_steer = mpc_vy + vy_aoa_err;

figure(39)
plot( u_steer, v_steer, 'vm'); grid on
hold on
plot(mpc_vy, mpc_ux, 'ob'); grid on

%% MSE wrt SNR and wrt number of antenna (M)
% disp('Get Vd output voltage')
% use_test_case = 1;
% az_deg = 20;
% el_deg = 10;
% if use_test_case
%     SNR_dB = 20;  
%     SNR = 10 .^( SNR_dB / 20 );
%     
%     R_m = 5; %radial distance of obj
%     az_obj = deg2rad( az_deg );  %polar distance of obj
%     el_obj = deg2rad( el_deg );
%     [ cylX_m, cylY_m, cylZ_m ] = Razel2xyz( R_m, az_obj, el_obj );
%     %[ checkU, checkV ] = xyz2uv_NT( cylX_m, cylY_m, cylZ_m, 0, 0 );
%     fig = figure;
%     plot(az_obj,el_obj,'vm');grid on
%     print(fig,'MySavedPlot','-dpng')
%     [ u_steer_0, v_steer_0 ] = azel2uv( az_obj, el_obj );
% 
% %   [ u_steer_0, v_steer_0 ] = azel2uv(az_obj,el_obj);
%     u_steer_0 = cos( az_obj ) .* sin( el_obj );
%     v_steer_0 = sin( el_obj );
% 
%     distCyl2Rx_X = cylX_m - array_x;
%     distCyl2Rx_Y = cylY_m - array_y;
%     distCyl2Rx_Z = cylZ_m - array_z;
%     
%     distCyl2Rx = sqrt( distCyl2Rx_X .^2 + distCyl2Rx_Y .^2 + distCyl2Rx_Z .^2 );        %% meters
%     
%     S21_store_sim = zeros( freq_sam, NelemTot );
%    
%     for ax = 1 : freq_sam
%         lambda = c / f_Hz( ax );
%         theor_phase_rad = mod( 2 .* pi .* distCyl2Rx ./ lambda, 2 * pi );
%         %Noise
%         Vt_new0( ax, : ) = SNR * exp( -1j * theor_phase_rad( : ) ).' + ( 1 / sqrt( 2 ) ) * ( randn( 1, NelemTot ) + 1j * randn( 1, NelemTot ) );
%     end
%    %Above gives output signal voltage
%     else
%     Vt_new0 = conj( Vt_new0 );
% end
%% 
disp('Use voltage for as')
x_signal = zeros (Nx, Ny, length(t));
x_sum = zeros (1,length(t));
x_signalstack = zeros (NSAs, Nx, Ny, length(t));
phase= zeros (Nx, Ny);   
array_x =zeros(Nx,Ny); array_y = zeros (Nx,Ny);

SA_steer = @(x,y,u,v,k)(exp(-1j*k*(x.*u + y.*v)));
atemp= [1:1:Nx];   %initial index, it starts from 1 not 0

theta_degree = 45; %30;   %angle_of_arrival
v_UAV = 100;       %UAV velocity in m/s
 f_Hz = 28e9; %2.4e9; 1e9

disp('Output Voltage Covariance') 
% For 1:NSAs, computeO/p of SAs for steerEachSA with phase comp for targ from each element of SA
for ai= 1: NSAs    

    array_x = arraystack_x (atemp, 1:Nx); %use atemp
    array_y = arraystack_y (atemp, 1:Ny); 
    
    steering_phase = SA_steer( array_x, array_y, -u_steer, -v_steer, k );
    
    atemp  = Nx*ai + [1:1:Nx]; %next loop value
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
    i
    j
    %complex valued analog signal using the I/Q sg at the pulse compressor input
   %steering_vector.Direction = [angle_of_arrival; 0];
    env_norm = @(lfm_scal, phi_init_rad, t, mu, wc, SNR)...
        (SNR * (exp ( 1j *(wc * t + mu * (t .^ 2 ) / 2 + phi_init_rad )) ) / lfm_scal ); 
    comp_env_norm = env_norm (lfm_scal,...
      phi_init_rad, t, mu, wc, SNR );
   
%% DOppler Shift
  
    doppler_shift = v_UAV*sin(theta_degree)*f_Hz/ c;
    doppler_beta  = exp(1i * 2 * pi * doppler_shift * (i+j));
    lfm_comp_env_norm = ...
    doppler_beta*comp_env_norm ;
      
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
    SAoutput_V( ai, : ) = tempSA; % NSAs x length(t) (in paper M x length(t))
end

Rv_est = (1/length(t)) * SAoutput_V * SAoutput_V'; %SAvoltage cov
figure(105)
plot(chirp)
grid on 
hold on
plot(SAoutput_V(:,1).')
plot(comp_env_norm)
%% Plot real(chirp) duration UAV speed
% Display results
%disp(['True UAV Speed: ' num2str(true_speed) ' m/s']);
%disp(['Estimated UAV Speed: ' num2str(estimated_speed) ' m/s']);

% Plot results
figure;
%subplot(2, 1, 1);
plot(real(chirp));
%plot(0:1/fs:size(chirp), real(chirp));

title('Transmitted Chirp Signal');
xlabel('Time (s)');
ylabel('Amplitude');
%%%
figure;
%subplot(2, 1, 1);
plot(imag(chirp));
%plot(0:1/fs:size(chirp), real(chirp));
title('Transmitted Chirp Signal');
xlabel('Time (s)');
ylabel('Imaginary Part');

%%
%% SAvoltages for the center of all the SAs. Next, BFout for this SA
display('Steer the array')
%mainbeam_3dB    =lambda/(Nelem*lambda);
%mainbeam_3dB_0  =0.0015;%0.02%0.1

mainbeam_3dB_0 = mainbeam_3dB / ( Bfrac / 8 );
x_3dB_0 = [ mainbeam_3dB_0/4 0 -mainbeam_3dB_0/4 0 ].';
y_3dB_0 = [ 0 mainbeam_3dB_0/4 0 -mainbeam_3dB_0/4 ].';
Nbeams = length(x_3dB_0);

Nbeams
beam_steer_ux = bsxfun( @plus, u_steer, x_3dB_0);
beam_steer_vy = bsxfun( @plus, v_steer, y_3dB_0);
%%
% display('Compute B matrix')
% SAstack = zeros (NSAs, Nbeams ) ;
% 
% for bx = 1:Nbeams
% atemp2 = [1:Nx];   %initial
% for ax= 1: NSAs % vary the SA from 1 to NSAs
% tempx   = arraystack_x (atemp2, 1:Nx); %use atemp
% tempy   = arraystack_y (atemp2, 1:Ny); 
% SAstack(ax, bx) =...
% (1/sqrt(Nelem_SA))...
% *sum( SA_steer( tempx(:), tempy(:), ...
% u_steer - beam_steer_ux( bx ), ...
% v_steer - beam_steer_vy( bx ), k) );
%         
% atemp2  = Nx*ax + [1:Nx]; %next loop value   %Nx+atemp2;
% end
% end

%Rinv= inv(Rv_est); %inverse of SAvoltage cov

%% Next is 2nd stage
%% specific funBased beam:fun is the search in spec dixn
display('Create beam and find its rank')
display('Compute B matrix')
Nbeams= 4; %8

SA_vec = zeros( NSAs, Nbeams );
B = zeros( Nbeams, NSAs );
for ax = 1 : Nbeams
    for bx = 1 : NSAs
        tempX = squeeze( SA_x_m( bx, :, : ) );
        tempY = squeeze( SA_y_m( bx, :, : ) );
        SA_vec( bx, ax ) = ( 1 / sqrt( NelemSA ) ) * sum( SA_steer( tempX( : ), tempY( : ), beam_steer_ux( ax ) - u_steer, beam_steer_vy( ax ) - v_steer, k ) );
        %SA_vec( bx, ax ) = ( 1 / sqrt( NelemSA ) ) * sum( SA_steer( tempX( : ), tempY( : ), beam_steer_ux( ax ), beam_steer_vy( ax ), k ) );
    end
end

for ax = 1 : Nbeams
    tempS = squeeze( SA_vec( 1 : NSAs, ax ) );
    B( ax, : ) = tempS / ( tempS' * tempS );
end


Rsam = SA_data_steered_out_vec.' * conj( SA_data_steered_out_vec );
% check line396 to line 404=ok 
for iBeam = 1 : Nbeams
    for iSA = 1 : NSAs
tempSA = SA_steer( arraystack_x( ( iSA - 1 ) * Nx + 1 : iSA * Nx, 1 : Ny ), arraystack_y( ( iSA - 1 ) * Ny + 1 : iSA * Ny, 1 : Ny ), u_steer - beam_steer_ux( iBeam ), ...
        v_steer - beam_steer_vy( iBeam ), k );
        Sn( iSA, 1 ) = ( 1 / ( Nx * Ny ) ) * sum( tempSA( : ) );
end
Bn( iBeam, 1 : NSAs ) = ( ( Rinv * Sn ) / ( Sn' * Rinv * Sn ) )';
end
% disp( ['Rank of estimated beamformer matrix = ', num2str( rank( Bn ) ) ] );

%% Employ Beam Bn and Rvest
display('For one time stamp,Compute cost function Employ Beam Bn and Rvest')
Snbr = 815; %one time stamp 
%covariance and its inv
Cn   = Bn * Rv_est * Bn';%beamVoltage cov
nu_2 = 0;
Cn_inv = inv( Cn + nu_2 * eye( iBeam ) );

%snapshot voltage in o/p beam
psi_0 = reshape( SAoutput_V( 1 : NSAs, Snbr ), NSAs, 1 );

%SA_cost:cost = 1 / real( SIR )
SA_cost_function_1 = ...
@( x ) SA_cost_function( x, k, ...
                       NSAs, Bn, Cn_inv, arraystack_x, arraystack_y, ...
                       psi_0, u_steer, v_steer, Nx, Ny );

x0 = [ u_steer, v_steer ];

%estimate the target AoA phase(ML) and fval/costout @ML 
options = optimset( 'MaxIter', 2e4, 'Display', 'iter', 'GradObj', 'off', 'Hessian', 'off', 'DerivativeCheck', 'off', 'TolFun', 1e-32, 'TolX', 1e-32, 'MaxFunEvals', 2e5 );

% [ ML, psi_out, exitflag ] = fminsearch( SA_cost_function_1, x0, options );
% % [ ML, psi_out, exitflag ] = fminunc( SA_cost_function_1, x0, options );

[ ML, cost_out, exitflag ] = fminsearch( SA_cost_function_1, x0, options );

mpc_ux
mpc_vy
ML
%% The N ×1 est beam vector a(un,vn) can be expressed as,
%a(un,vn)=(1/?M)*e?j*(2?/lambda)*(c_(x,m)d_x*u_n+c_(y,m)d_y*v_n), (0 ? n ? N ? 1).
u_est = ML( 1 ) 
v_est = ML( 2 ) 
u_steer
v_steer

u_steer_0
v_steer_0

u_error = u_est - u_steer_0
v_error = v_est - v_steer_0

u_error_0 = u_steer - u_steer_0
v_error_0 = v_steer - v_steer_0

rms_error = sqrt( u_error^2 + v_error^2 )

rms_error_0 = sqrt( u_error_0^2 + v_error_0^2 )


figure(109)
plot(beam_steer_ux,beam_steer_vy,'+g')  %Rosette Beams
grid on
hold on
plot(mpc_ux,mpc_vy,'or')                %Targt = TRUTH
plot(u_steer,v_steer,'+b')              %Steer = 
plot(ML(1),ML(2),'ok')
hold off
legend('Rosette Beams','Truth','Steer','Estimate');
toc
%%
figure(212)  % subplot(2, 1, 2);
plot(1:num_elements_x, estimated_speed); %wo BF
plot(1:Nx, ML);
 
title('Estimated UAV Speed for Each Array Element');
xlabel('Element Index');
ylabel('Estimated Speed (m/s)');
chirp_signal_baseband = chirp(0:1/fs:duration, baseband_frequency, duration, chirp_bandwidth);

if show_all_plots

BandWidth_Hz = 0 : 1e9 : 14e9;
Nf = length( BandWidth_Hz );

for ax = 1 : Nf
    [ rms_error( ax ), rms_error_0( ax ) ] ...
        = fast_beamspace_MLE_version( BandWidth_Hz( ax ) );
end
%  average_error = mean( [rms_error] ) 
average_error = mean( [0.01 0.017 0.028 0.042 0.012 0.022 0.028 0.034 0.05 ...
                       0.018 0.06 0.059 0.022 0.024 0.028] )
                   
%%
% rms_error     rms_error_0
% 'New Error', 'Original Error'

figure( 12 )
plot( BandWidth_Hz / 1e9, rms_error, '-o', 'LineWidth', 1.6, 'MarkerSize', 6, 'MarkerFace', [ 0.2 0.36 0.99 ] );
hold on
plot( BandWidth_Hz / 1e9, rms_error_0, '-o', 'LineWidth', 1.6, 'MarkerSize', 6, 'MarkerFace', [ 0.02 0.36 0.79 ] );
hold off
grid on
xlabel('Bandwidth [GHz]');
ylabel('RMS AOA Error');
eval( [ 'title( ''RMS AOA Error vs Bandwidth'', ''FontSize'', 14, ''Interpreter'', ''Latex'' )' ] );
legend( 'New Error', 'Original Error', 'Location', 'South' );
if print_plots
    set( gcf, 'PaperPositionMode', 'auto' );
    eval( [ 'print( ''-dpng'', [ PlotPath, ''error_vs_bandwidth.png''] )' ] );
    pause( 3 )
end

end
%%
if show_all_plots_2
    
%Snbr = 311;
[ I1, I2 ] = max( 20 * log10( abs( SA_data_steered_out_vec( : ) ) ) );
[ Bin_nbr, SA_nbr ] = ind2sub( [ freq_sam NSAs ], I2 );
xi = B * SA_data_steered_out_vec.';
psi_0 = xi( :, Bin_nbr );


check_beam = compute_FFT( UplotFFT1, VplotFFT1, squeeze( SA_vec( :, 4 ) ), phase_center_X, phase_center_Y, lambda );
check_beam_dB = 20 * log10( abs( check_beam ) );
figure( 462 )
pcolor( UplotFFT1, VplotFFT1, check_beam_dB );
xlabel( 'U' )
ylabel( 'V' );
title('Rosette Beam');
colormap jet; shading flat; axis('square'); colorbar; %caxis([ D11-65 D11 ])
if print_plots
    set( gcf, 'PaperPositionMode', 'auto' );
    eval( [ 'print( ''-dpng'', [ PlotPath, ''rosette_beam.png''] )' ] );
    pause( 1 );
end


%% Composite Beam
display('Start Beam of Beams');

%C = B * Rsam * B';
C = B * B';

rank_C = rank( C )
rank_B = rank( B )

%C_inv = eye( Nbeams );
nu_3 = 0;
C_inv = inv( C + nu_3 * eye( Nbeams ) );

%%
Duv = 0.5;  %0.1  %0.5;     %0.05;
Umin = -Duv;  
Vmin = -Duv; 
Umax = Duv;  
Vmax = Duv;  

Ndelta = 100;  %200; 
cost_func = zeros( Ndelta, Ndelta );
ugrid = linspace( Umin, Umax, Ndelta );
vgrid = linspace( Vmin, Vmax, Ndelta );

[ DeltaUgrid, DeltaVgrid ] = meshgrid(  ugrid, vgrid  );


%% cost_function
display('Compute cost function');

for jx = 1 : Ndelta2
    for lx = 1 : Ndelta2
        uv00 = [ u_steer + DeltaUgrid2( jx, lx ), v_steer + DeltaVgrid2( jx, lx ) ];
        cost_func2( jx, lx ) = SA_cost_function_data( uv00, k, NSAs, B, C_inv, SA_x_m, SA_y_m, psi_0, u_steer, v_steer );
    end
    jx
end
cost_func_dB2 = 20 * log10( abs( cost_func2 ) );



%%
figure( 122 )
plot3( u_steer + DeltaUgrid, v_steer + DeltaVgrid, cost_func_dB )
grid on
xlabel( 'U' ); ylabel( 'V' ); zlabel( 'Cost (dB)' );
title( 'Pattern of Composite Beam' )
view( 62, 80 )
if print_plots
    set( gcf, 'PaperPositionMode', 'auto' );
    eval( [ 'print( ''-dpng'', [ PlotPath, ''cost_mle.png''] )' ] );
    pause( 1 )
end


L1 = min( cost_func_dB( : ) ) - 5;
L2 = max( cost_func_dB( : ) ) + 5;
figure( 123 )
pcolor( u_steer + DeltaUgrid, v_steer + DeltaVgrid, cost_func_dB )
grid on
xlabel( 'U' ); ylabel( 'V' ); zlabel( 'Cost (dB)' );
title( 'Pattern of Composite Beam' )
colormap jet; colorbar; caxis([ L1 L2 ]);  %caxis([ min( cost_func_dB( : ) ) max( cost_func_dB( : ) ) ])
if print_plots
    set( gcf, 'PaperPositionMode', 'auto' );
    eval( [ 'print( ''-dpng'', [ PlotPath, ''cost_mle_pcolor.png''] )' ] );
    pause( 1 )
end

%%
figure( 124 )
plot3( u_steer + DeltaUgrid2, v_steer + DeltaVgrid2, cost_func_dB2 )
grid on
xlabel( 'U' ); ylabel( 'V' ); zlabel( 'Cost (dB)' );
title( 'Pattern of Composite Beam' )
view( 62, 80 )
if print_plots
    set( gcf, 'PaperPositionMode', 'auto' );
    eval( [ 'print( ''-dpng'', [ PlotPath, ''cost_mle_v2.png''] )' ] );
    pause( 1 )
end


L11 = min( cost_func_dB2( : ) ) - 5;
L22 = max( cost_func_dB2( : ) ) + 5;
figure( 125 )
pcolor( u_steer + DeltaUgrid2, v_steer + DeltaVgrid2, cost_func_dB2 )
grid on
xlabel( 'U' ); ylabel( 'V' ); zlabel( 'Cost (dB)' );
title( 'Pattern of Composite Beam' )
colormap jet; colorbar; caxis([ L11 L22 ]);  %caxis([ min( cost_func_dB( : ) ) max( cost_func_dB( : ) ) ])
if print_plots
    set( gcf, 'PaperPositionMode', 'auto' );
    eval( [ 'print( ''-dpng'', [ PlotPath, ''cost_mle_pcolor_v2.png''] )' ] );
    pause( 1 )
end
% figure( 126 )
% plot( beam_steer_ux, beam_steer_vy, 'o', 'MarkerSize', 7, 'MarkerFace', [ 0.02 0.36 0.49 ] );
% hold on
% plot( u_steer, v_steer, 'rx', 'MarkerSize', 7 )
% plot( u_steer_0, v_steer_0, 'kx', 'MarkerSize', 7 );
% plot( u_est, v_est, 'bx', 'MarkerSize', 7 );
% hold off
% grid on
% xlabel('U');
% ylabel('V');
% eval( [ 'title( ''Rosette Beam Locations:  Bandwidth = ', num2str( BandWidth_Hz / 1e9 ), ' GHz'', ''FontSize'', 14, ''Interpreter'', ''Latex'' )' ] );
% legend('Rosette', 'Steer', 'Object', 'Estimate' );
% if print_plots
%     set( gcf, 'PaperPositionMode', 'auto' );
%     eval( [ 'print( ''-dpng'', [ PlotPath, ''rosette_beams_', num2str( BandWidth_Hz / 1e9 ), '.png''] )' ] );
%     pause( 3 )
% end

end






















close all
clear all
clc
%%
tic
fs = 8e9;          %  sampling frequency
delF = 1e9;
fc = 0;             %  baseband frequency
wc = 2 * pi * fc;
Ts = 1 / fs;        %  sampling period
T = 1e-6;  %10e-6;           %  chirp duration
t = 0 : Ts : T;
Nsam = length( t );

SNR_dB = 60;
phi_init_rad = 0;
BT = delF * T;          %  time-bandwidth product
SNR = 10 ^ ( SNR_dB / 20 ); %snr in voltage ratio
lfm_scal = 1;
%  Calculate frequency sweep rate (radians/sec)
mu = 2 * pi * delF / T;
%% Beam of beam of subarrays 1. the kth element of array belongs to the sth subarray the ith lement
Nx0= 30;
Ny0= 30;
Nelem= Nx0 * Ny0;
NSAs = 9;
% Assuming 1st SA at center (0,0) (position in matrix form is 1st= 11)
Nx= 10; %sqrt (Nelem /NSAs);
Ny = 10;  %sqrt (Nelem/ NSAs);
Nelem_SA = Nx * Ny;
%%
% f_Hz =1e9;
f_Hz =28e9;  
c=3e8; lambda = c/f_Hz; k = 2*pi/lambda; spacing = lambda/2; % warning('off','all');

%% clear tempx      %square NSAs_x =NSAs_y ; Nx = Ny
NSAs_x = 3;     NSAs_y = NSAs/ NSAs_x;
a_x = zeros(NSAs_x* Nx, Nx);
a_y = zeros(NSAs_y* Ny, Ny);

temp = [1:Nx];           %initial index, it starts from 1 not 0
g = [0:Nx-1];            % initial elemt position in a SA, starts from 0

for ia_x = 1: NSAs_x    % combine SA matx from 1 to NSAs_x
% a_x( temp, 1:Nx ) = kron((grid),ones (Nx,1)) ;      %only a_x

[a_x(temp, 1:Nx) a_y(temp, 1:Ny)] = meshgrid(g , g) ;

temp  = Nx + temp  ;    %next loop value %   temp_x = Nx * ia_x + [1:Nx];
    g = Nx + g;
end

array_x = a_x * spacing; %30x10 %1st 2nd 3rd x/y of SA, rwise
array_y = a_y * spacing ;

%% stack o/ps % stack periodic matrix row wise   % perm = repmat(eye(NSAs_x), [NSAs_x 1]);

perm = repmat( eye(NSAs_x * Nx), [NSAs_x 1] );
stack_a_x = perm * a_x;         %90x10 %all SA, rwise
arraystack_x = stack_a_x * spacing;

stack_a_y = perm * a_y;
arraystack_y = stack_a_y * spacing;

%% plot all the SA in x-y plane
 plot_array = 1;
%plot_array = 0;
if plot_array
for short = 1:1
    
origin_vec =  [0: Nx-1];
origin_pos =  [-1 0 0; 0 0 0; 0 0 0];
origin_next_pos = eye (NSAs_x) +  origin_pos;

initial_index = Nx * origin_next_pos * [0: NSAs_x-1 ].';
next_increment= repmat( initial_index, [1, Nx]  );

order_xy = repmat(origin_vec, [NSAs_x, 1])+ next_increment;
order_SA = order_xy* spacing;

[xall_SA yall_SA] = meshgrid( order_SA, order_SA);

%% plot max min cent0 of SA

atemp = [1:Nx];         %initial index, it starts from 1 not 0
ia = NSAs ;
SA_x0 =zeros (NSAs,1) ;  
SA_y0 = zeros(NSAs,1);  array_x0= zeros(NSAs) ; array_y0 =zeros (NSAs);
%%
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
%%
gtemp = [1:Nx];
array_x11 = arraystack_x( gtemp, gtemp  );
array_x12 = arraystack_x( Nx + gtemp, gtemp  );
array_x13 = arraystack_x( 2*Nx + gtemp, gtemp );

array_y11 = arraystack_y( gtemp, gtemp  );
array_y12 = arraystack_y( gtemp, gtemp );
array_y13 = arraystack_y( gtemp, gtemp );

array_x21 = arraystack_x( gtemp, gtemp ); 
array_x22 = arraystack_x( Nx + gtemp, gtemp );
array_x23 = arraystack_x( 2*Nx + gtemp, gtemp );

array_y21 = arraystack_y( Nx + gtemp, gtemp );
array_y22 = arraystack_y( Nx + gtemp, gtemp );
array_y23 = arraystack_y( Nx + gtemp, gtemp );

array_x31 = arraystack_x( gtemp, gtemp );
array_x32 = arraystack_x( Nx + gtemp, gtemp );
array_x33 = arraystack_x( 2*Nx + gtemp, gtemp );

array_y31 = arraystack_y( 2*Nx + gtemp, gtemp );
array_y32 = arraystack_y( 2*Nx + gtemp, gtemp );
array_y33 = arraystack_y( 2*Nx + gtemp, gtemp );

   show_array = 1;         %change condition to show or not
    %%
    if show_array
    figure(99)
    plot( SA_x00, SA_y00, 'vk'); grid on
    hold on
    plot(array_x11,array_y11,'ob');
    plot(array_x12,array_y12,'ok');
    plot(array_x13,array_y13,'om');
    plot(array_x21,array_y21,'oc');
    plot(array_x22,array_y22,'og');
    plot(array_x23,array_y23,'or');
    plot(array_x31,array_y31,'ob');
    plot(array_x32,array_y32,'oy');
    plot(array_x33,array_y33,'om');
    hold off
    grid on
    end
end
end
%% phase computation of target location& add possible Error for each element
L = length(t);
rng(0);

targY = 7e3/3.28084;          %[7e3 10.5e3 13e3 9.5e3 14.5e3]/3.28084;
targX = 3 * -14e3/3.28084;    %3 * [-14e3 -7e3 2e3 7e3 13e3]/3.28084;
targR = 200 * 1e2;            %[200 180 197.5 175 199.5] * 1e2;

%find locasn f trgt z-dim (m) from given range, x, and y (objectLocation)
targZ = sqrt(targR.^2 - targX.^2 - targY.^2);
targPhi = atan2(targY,targX);           % target angle phi from given x and y
targTh = acos(targZ./targR);            % target angle theta from given z and range
[ mpc_ux, mpc_vy ] = fThetaphi2uv( targTh, targPhi );

%%
%UserDefinedFormulaError bw steeringDixn n trueSigDixn in 1dim(x) 
Bfrac = 10;
mainbeam_3dB = 0.886 * lambda / ( Nx0 * lambda );   

ux_aoa_err = ( mainbeam_3dB / Bfrac );      %rand( 1, array.nTarg );        %  user defined error between steering direction and true signal direction %( mainbeam_3dB / 8 )*[1 -1 1-1] in quadrant beam error
vy_aoa_err = ( mainbeam_3dB / Bfrac );      %rand( 1, array.nTarg );
initial_rms_error = sqrt( ux_aoa_err .^2 + vy_aoa_err .^2 );
%error

u_steer = mpc_ux + ux_aoa_err;            %recv beam steering direction
v_steer = mpc_vy + vy_aoa_err;

%%
x_signal = zeros (Nx, Ny, length(t));
x_sum = zeros (1,length(t));
x_signalstack = zeros (NSAs, Nx, Ny, length(t));
phase= zeros (Nx, Ny);   
array_x =zeros(Nx,Ny); array_y = zeros (Nx,Ny);

SA_steer = @(x,y,u,v,k)(exp(-1j*k*(x.*u + y.*v)));
atemp= [1:1:Nx];   %initial index, it starts from 1 not 0

% For 1:NSAs, computeO/p of SAs for steerEachSA with phase comp for targ from each element of SA
for ai= 1: NSAs    

    array_x = arraystack_x (atemp, 1:Nx); %use atemp
    array_y = arraystack_y (atemp, 1:Ny); 
    
    steering_phase = SA_steer( array_x, array_y, -u_steer, -v_steer, k );
    
    atemp  = 10*ai + [1:1:Nx]; %next loop value
    tempSA = zeros (length(t), 1);

    for i= 1 : Nx  % Nx=Ny=10 %1st SA
    for j= 1 : Ny
    %[xQ yQ zQ] = fAzel2xy(azQ, elQ, r); %1.object position
    xQ = targX;
    yQ = targY;
    zQ = targZ;
   
    dx =xQ - array_x (i,j) ;%2. distance bw object and every element
    dy =yQ - array_y (i,j);
    dz =zQ;
    d = sqrt(dx.^2 + dy.^2 + dz.^2 ); %

    phase(i,j) = mod( (2*pi/lambda)* d, 2 * pi );%3.1. phase_shift based on pathshift
    phi_init_rad = phase(i,j);%3.2 chirp for this phase ip

    %complex valued analog signal using the I/Q sg at the pulse compressor input
   
  env_norm = @(lfm_scal, phi_init_rad, t, mu, wc, SNR)...
        (SNR * (exp ( 1j *(wc * t + mu * (t .^ 2 ) / 2 + phi_init_rad )) ) / lfm_scal ); 
  lfm_comp_env_norm = env_norm (lfm_scal, phi_init_rad, t, mu, wc, SNR );
    
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
%% SAvoltages for the center of all the SAs. Next, BFout for this SA

mainbeam_3dB_0 = mainbeam_3dB / ( Bfrac / 8 );
x_3dB_0 = [ mainbeam_3dB_0/4 0 -mainbeam_3dB_0/4 0 ].';
y_3dB_0 = [ 0 mainbeam_3dB_0/4 0 -mainbeam_3dB_0/4 ].';
Nbeams = length(x_3dB_0);

beam_steer_ux = bsxfun( @plus, u_steer, x_3dB_0);
beam_steer_vy = bsxfun( @plus, v_steer, y_3dB_0);
%%
SAstack = zeros (NSAs, Nbeams ) ;

for bx = 1:Nbeams
    atemp2 = [1:Nx];   %initial
    for ax= 1: NSAs % vary the SA from 1 to NSAs
        tempx = arraystack_x (atemp2, 1:Nx); %use atemp
        tempy = arraystack_y (atemp2, 1:Ny); 
        SAstack(ax, bx) = (1/sqrt(Nelem_SA))...
            *sum( SA_steer( tempx(:), tempy(:), ...
            u_steer - beam_steer_ux( bx ), v_steer - beam_steer_vy( bx ), k) );
        
        atemp2  = 10*ax + [1:Nx]; %next loop value   %Nx+atemp2;
    end

end

Rinv= inv(Rv_est); %inverse of SAvoltage cov

%% matrix of weight vectors
B = zeros(Nbeams, NSAs); % 4x9
for bx = 1:Nbeams
    w = Rinv * SAstack(:,bx) /(SAstack(:, bx)' * Rinv * SAstack(:, bx));
    B(bx,:) = w';

end
%% xi,C are op of 1st stage 
xi =B * SAoutput_V;
C = B * Rv_est * B';

%% Next is 2nd stage
%% specific funBased beam:fun is the search in spec dixn
Nbeams= 4;
for iBeam = 1 : Nbeams
    for iSA = 1 : NSAs
        tempSA = SA_steer( arraystack_x( ( iSA - 1 ) * Nx + 1 : iSA * Nx, 1 : Ny ), arraystack_y( ( iSA - 1 ) * Ny + 1 : iSA * Ny, 1 : Ny ), u_steer - beam_steer_ux( iBeam ), ...
            v_steer - beam_steer_vy( iBeam ), k );
        Sn( iSA, 1 ) = ( 1 / ( Nx * Ny ) ) * sum( tempSA( : ) );
    end
    Bn( iBeam, 1 : NSAs ) = ( ( Rinv * Sn ) / ( Sn' * Rinv * Sn ) )';
end
% disp( ['Rank of estimated beamformer matrix = ', num2str( rank( Bn ) ) ] );

%%
Snbr = 713; %one time stamp 

Cn = Bn * Rv_est * Bn';%beamVoltage cov
nu_2 = 0;
Cn_inv = inv( Cn + nu_2 * eye( iBeam ) );

%snapshot voltage in o/p beam
psi_0 = reshape( SAoutput_V( 1 : NSAs, Snbr ), NSAs, 1 );

%SA_cost:cost = 1 / real( SIR )
SA_cost_function_1 = @( x ) SA_cost_function( x, k, ...
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

figure(103)
plot(beam_steer_ux,beam_steer_vy,'+g')  %Rosette Beams
grid on
hold on
plot(mpc_ux,mpc_vy,'or')                %Targt = TRUTH
plot(u_steer,v_steer,'+b')              %Steer = 
plot(ML(1),ML(2),'ok')
hold off
legend('Rosette Beams','Truth','Steer','Estimate');
toc

%%
tic
clear all
close all
PlotPath = 'C:\Users\Abhishek Agrahari\Documents\MVDR\Plots_vs_Elements_060224\Channel_Doppler_Shift\UAV\fig_estimation_UAV\'; 
show_all_plots =0;
show_all_plots_2 =1;

vd = 3e1;
theta = pi/6;
freq_sam =1351;
Bandwidth_Hz =78e6;
fLO = 2.4e9;
fHI = fLO + Bandwidth_Hz;
fMID = 0.5*(fLO+fHI);
%Input range
f_Hz      = linspace( fLO, fHI, freq_sam ); 
c         = 3e8;
lambdaHI  = c / fHI; 
spacing   = lambdaHI / 2;
lambdaMID = c / fMID;  k = 2 * pi / lambdaHI;
lambda    = c ./ f_Hz;
%
fs        = 8e9;          %sampling frequency
delF      = 1e9;
fc        = 2.4e9;         %baseband frequency
wc        = 2 * pi * fc;
Ts        = 1 / fs; %  sampling period
T         = 1e-6;  %10e-6;           %  chirp duration
t         = 0 : Ts : T;
Nsam      = length( t ); % observation time
%%
SNR_dB      = 2 %0.02;%5;20;60;
phi_init_rad= 0;
BT          = delF * T;          %  time-bandwidth product
SNR         = 10 ^ ( SNR_dB / 20 ); %snr in voltage ratio
lfm_scal    = 1;
%%  Calculate frequency sweep rate (radians/sec)
mu          = 2 * pi * delF / T;
Array_steer = @(x,y,u,v,k)( exp(-1j*k*(x.*u + y.*v)) );
%% Phased array GS Elements
Nx0         = 2;
Ny0         = 2;
Nelem       = Nx0*Ny0; %36
 
NelemOver2  = ( Nelem - 1 ) / 2; 
NelemTot    = Nelem ;
%NelemTot   = Nelem*Nelem ;
NSAs        = 1;
save_array  = 1;
steer_off   = 1;
%%
use_with_FFT =1;
Nbeams =4;

if use_with_FFT
    disp('array layout')
[x_element_index, y_element_index ]...
    =meshgrid(0:Nx0-1, 0:Ny0-1);
array_x_planar...
=x_element_index *spacing;
array_y_planar...
=y_element_index*spacing;
array_z_planar...
    =zeros(size(array_x_planar));
x_element_index = reshape( x_element_index, 1, NelemTot );
y_element_index = reshape( y_element_index, 1, NelemTot );
rng(0);
use_obj_case = 1;
if use_obj_case
    R_m =5;
    az_obj=deg2rad(20);
    el_obj=deg2rad(10);
    
    [cylX_m, cylY_m, cylZ_m] = Razel2xyz( R_m, az_obj, el_obj );
    distCyl2Rx_X = cylX_m - array_x_planar;
    distCyl2Rx_Y = cylY_m - array_y_planar;
    distCyl2Rx_Z = cylZ_m - array_z_planar;
    distCyl2Rx   = sqrt( distCyl2Rx_X .^2 + distCyl2Rx_Y .^2 + distCyl2Rx_Z .^2 );        %% meters
    [checkU,checkV] = xyz2uv_NT(cylX_m, cylY_m, cylZ_m,0,0);
    
    [u_steer_0, v_steer_0]...
    = azel2uv(az_obj,el_obj);
   S21_store_sim = ...
       zeros(freq_sam, NelemTot);    

   for ax  = 1 : freq_sam
    lambda = c / f_Hz( ax );  
    theor_phase_rad  = mod( 2 .* pi .* distCyl2Rx ./ lambda, 2 * pi );
    %Noise
    Vt_new0( ax, : ) = SNR * ...
        exp( -1j * theor_phase_rad( : ) ).' ...
        + ( 1 / sqrt( 2 ) ) * ( randn( 1, NelemTot )...
        + 1j * randn( 1, NelemTot ) );
   end
    else
    Vt_new0         = conj( Vt_new0 );
end
Vt_new_rank         = rank( Vt_new0 );
%Vt_new               = ifft( Vt_new0 );
    Nbeams              = 4 %16;  %8;
    Nbeams
%%
if steer_off
    v_steer = v_steer_0 -  0.0152;
    u_steer = u_steer_0 - -0.0142;
else
    v_steer = v_steer_0;
    u_steer = u_steer_0;
end
%%
indicator = 1
if (indicator==1)
Doppler_shift=(vd*sin(theta)*fc)/(c)
else
Doppler_shift=-(vd*sin(theta)*fc)/(c)
end
fd    = Doppler_shift;
beta1 = (2*pi*fd).*[0:1:(Nbeams-1)];
complexj = -j*beta1*Ts;
Doppler_matrix_beta = diag([exp(complexj)].');
mainbeam_3dB   = lambda /(Nelem*lambda/2);
mainbeam_3dB_0 = mainbeam_3dB;
      %end
%end
%%
DeltaBeam   = 2*pi/Nbeams;

BeamSpacing = linspace(0,2*pi-DeltaBeam, Nbeams);
beamU       = cos( BeamSpacing );
beamV       = sin( BeamSpacing );
Doppler    = Doppler_shift*Ts*linspace(0,Nbeams-1,Nbeams);

shift_beamU = 0.5*Doppler.'+mainbeam_3dB_0*beamU.';
shift_beamV = 0.5*Doppler.'+mainbeam_3dB_0*beamV.';

%% Add phase change due to Dopppler
beam_steer_ux1...
    = bsxfun(@plus,u_steer, shift_beamU);
beam_steer_vy1...
    = bsxfun(@plus,v_steer,shift_beamV);
disp('shifted steer beam and wo shifted')
wo_Doppler_ux1...
=bsxfun(@plus, u_steer, mainbeam_3dB_0*beamU.');
wo_Doppler_vy1...
=bsxfun(@plus, u_steer, mainbeam_3dB_0*beamV.');
wo_steer_ux1...
=mainbeam_3dB_0*beamU.';
wo_steer_vy1...
=mainbeam_3dB_0*beamV.';
%%
disp('Process EntireArray')
steering_phase...
= Array_steer( array_x_planar,...
array_y_planar, ...
u_steer, ...
v_steer, k );
steered_array_data...
= bsxfun( @times,...
steering_phase( : ).', ...
Vt_new0 );
steered_array_data...
=  Vt_new0;
%%
disp('FFT of steered array data ')
Nfft = 256;
DeltaU_Nfft...
=[-Nfft/2:Nfft/2-1]*(2/Nfft);
DeltaV_Nfft...
    =[-Nfft/2:Nfft/2-1]*(2/Nfft);
UplotFFT = zeros(Nfft,Nfft);
VplotFFT = zeros(Nfft,Nfft);
for lx = 1:Nfft
    for nx = 1:Nfft
UplotFFT( lx, nx ) = DeltaU_Nfft( lx );
VplotFFT( lx, nx ) = DeltaV_Nfft( nx );    
    end
end
R2 =sqrt(UplotFFT.*UplotFFT+VplotFFT.*VplotFFT);

UplotFFT1 = UplotFFT;
VplotFFT1 = VplotFFT;

UplotFFT1(R2>1)=NaN;
VplotFFT(R2>1) =NaN;

check_FTout = compute_FFT(...
    UplotFFT1, VplotFFT1,...
    squeeze( steered_array_data( ...
    freq_sam, : ) ).', array_x_planar( : ), array_y_planar( : ), lambda );
check_FTout_dB = 20*log10(abs( check_FTout ));
%%
print_plots = 1
display('steered array plots for fix position')
figure( 412 )
pcolor( UplotFFT1, VplotFFT1, check_FTout_dB );
xlabel( 'U' )
ylabel( 'V' );
title('Steered Array Output');
colormap jet; shading flat; axis('square'); colorbar; %caxis([ D11-65 D11 ])
%         if print_plots
%             set( gcf, 'PaperPositionMode', 'auto' );
%             eval( [ 'print( ''-dpng'', [ PlotPath, ''steered_array_output.png''] )' ] );
%             pause( 1 );
%         end
%%
disp('Compute B Matrix')

check_FTout = compute_FFT( UplotFFT1, VplotFFT1, squeeze( steered_array_data0( freq_sam, : ) ).', array_x_planar( : ), array_y_planar( : ), lambda );
check_FTout_dB = 20 * log10( abs( check_FTout ) );

%%
display('Phase Center');
Mx=Nx0;
My=Ny0; 
for mx = 1 : Nelem
for nx = 1 : Nelem
Xmin = min( x_element_index );
Xmax = max( x_element_index );
Ymin = min( y_element_index );
Ymax = max( y_element_index );

[ SAx, SAy ] = meshgrid( -Mx/2 : Mx/2 - 1, -My/2 : My/2 - 1 );

x_ind_SA = SAx + mx + Mx/2 - 1;
y_ind_SA = SAy + nx + My/2 - 1;
                
if( ( max( max( x_ind_SA + Xmin ) ) <= Xmax ) && ( min( min( x_ind_SA + Xmin ) ) >= Xmin ) && ( min( min( y_ind_SA + Ymin ) ) >= Ymin ) && ( max( max( y_ind_SA + Ymin ) ) <= Ymax ) )
SA_index_X_store( cnt, 1 : Mx, 1 : My ) = x_ind_SA;
SA_index_Y_store( cnt, 1 : Mx, 1 : My ) = y_ind_SA;

SA_x_m( cnt, 1 : Mx, 1 : My ) = ( Xmin + x_ind_SA ) * spacing;
SA_y_m( cnt, 1 : Mx, 1 : My ) = ( Ymin + y_ind_SA ) * spacing;

cnt = cnt + 1;
end
end
end
for ax = 1 : NSAs
    phase_center_X( ax, 1 ) = mean( mean( SA_x_m( ax, :, : ) ) );
    phase_center_Y( ax, 1 ) = mean( mean( SA_y_m( ax, :, : ) ) );
end
beam_steer_ux = beam_steer_ux1;
beam_steer_vy = beam_steer_vy1;
beam_steer_ux2= mainbeam_3dB_0*beamU.';
beam_steer_vy2= mainbeam_3dB_0*beamV.';

Array_vec = zeros(Nelem,Nbeams);
B      = zeros( Nelem, NSAs );

for ax = 1 : Nbeams
tempS= ( 1 / sqrt( Nelem ) )...
    * Array_steer( beam_steer_ux( ax ), beam_steer_vy( ax ), phase_center_X, phase_center_Y, k );
 
Array_vec( :,ax )  = tempS ;
 %B( ax, : )     = tempS.' / ( tempS' * tempS );
 B( ax, : )        = tempS.';
end











































end

















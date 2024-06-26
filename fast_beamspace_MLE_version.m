%%  This is a fast version of the MLE beamspace algorithm

function [ rms_error, rms_error_0 ] = fast_beamspace_MLE_version( BandWidth_Hz, PlotPath )




% clc
% clear all
% close all
% BandWidth_Hz = 0;
% unix('cp -Rf *.m BACKUP');
% unix('cp -Rf *.m BACKUP_AGAIN');
 

%% PlotPath = '/home/pgv/Davis_Fante_Simulation/With_Position_Errors/Fast_Version/Plots_041620_CASE_4/';
% PlotPath = '\\bd5\5G_Share\Abhishek\E_band_feb12\temp\MVDR\Matlab_Code_for_Abhishek_100219\Samurai_Conference_Room\BoB_overlapping_array\fig_array_shift_SA\';

show_all_plots = 0;
show_all_plots_2 = 0;

print_plots = 0;

steer_off = 3;

apply_errors = 0;  % 0, 1,  2 or 3
applyBF_error = 0;  %  0 or 1
apply_position_errors = 0;  %1;
SA_phase_center_error = 0;

SA_steer = @(x,y,u,v,k)( exp( 1j * k * ( x.*u + y.*v ) ) );

c = 3e8;
freq_sam = 1351; Nelem = 35; NelemOver2 = ( Nelem - 1 ) / 2; NelemTot = Nelem * Nelem;

%BandWidth_Hz = 0;  %10e9;  %1e9;
fLO = 26.5e9;
fHI = fLO + BandWidth_Hz;
fMID = 0.5 * ( fLO + fHI );

f_Hz = linspace( fLO, fHI, freq_sam ); 
lambdaHI = c / fHI; spacing = lambdaHI / 2;
lambdaMID = c / fMID;  k = 2 * pi / lambdaHI;

[ x_element_index, y_element_index ] = meshgrid( 0 : Nelem - 1, 0 : Nelem - 1 );
array_x_planar = x_element_index * spacing;
array_y_planar = y_element_index * spacing;
array_z_planar = zeros( Nelem, Nelem );
x_element_index = reshape( x_element_index, 1, NelemTot );
y_element_index = reshape( y_element_index, 1, NelemTot );

rng( 0 );

use_test_case = 1;

if use_test_case
    SNR_dB = 20;  
    SNR = 10 .^( SNR_dB / 20 );
    
    R_m = 5;
    az_obj = deg2rad( 20 );  
    el_obj = deg2rad( 10 );
    
    clear Vt_new
    
    [ cylX_m, cylY_m, cylZ_m ] = Razel2xyz( R_m, az_obj, el_obj );
    [ checkU, checkV ] = xyz2uv_NT( cylX_m, cylY_m, cylZ_m, 0, 0 );

    [ u_steer_0, v_steer_0 ] = azel2uv( az_obj, el_obj );

    distCyl2Rx_X = cylX_m - array_x_planar;
    distCyl2Rx_Y = cylY_m - array_y_planar;
    distCyl2Rx_Z = cylZ_m - array_z_planar;
    
    distCyl2Rx = sqrt( distCyl2Rx_X .^2 + distCyl2Rx_Y .^2 + distCyl2Rx_Z .^2 );        %% meters
    
    S21_store_sim = zeros( freq_sam, NelemTot );
   
    for ax = 1 : freq_sam
        lambda = c / f_Hz( ax );
        theor_phase_rad = mod( 2 .* pi .* distCyl2Rx ./ lambda, 2 * pi );
        %Noise
        Vt_new0( ax, : ) = SNR * exp( -1j * theor_phase_rad( : ) ).' + ( 1 / sqrt( 2 ) ) * ( randn( 1, NelemTot ) + 1j * randn( 1, NelemTot ) );
    end
else
    Vt_new0 = conj( Vt_new0 );
end
Vt_new_rank = rank( Vt_new0 )

XmidArray = ( max( array_x_planar( : ) ) + min( array_x_planar( : ) ) ) / 2;
YmidArray = ( max( array_y_planar( : ) ) + min( array_y_planar( : ) ) ) / 2;
array_distance_from_center_m = sqrt( ( array_x_planar - XmidArray ).^2 + ( array_y_planar - YmidArray ).^2 );

MaxPh_err_rad_3 = pi/4;
MaxAmp_err_dB_3 = 1;  %3;
MaxAmp_err_3 = 10 ^ ( MaxAmp_err_dB_3 / 20 );
Ph_err_rad_3 = MaxPh_err_rad_3 * sin( ( 2 * pi / lambdaHI ) * ( array_distance_from_center_m / 5 ) );
Amp_err_3 = MaxAmp_err_3 * sin( ( 2 * pi / lambdaHI ) * ( array_distance_from_center_m / 10 ) );

if apply_position_errors == 1
    complexPosErrors = repmat( reshape( ( 1 + Amp_err_3 ) .* exp( 1j *  Ph_err_rad_3 ), 1, NelemTot ), freq_sam, 1 );
else
    complexPosErrors = ones( freq_sam, NelemTot );
end

Vt_new0 = Vt_new0 .* complexPosErrors;
Vt_new = ifft( Vt_new0 );


%%
Mx = 7;
My = 7;
NelemSA = Mx * My;

Xmin = min( x_element_index );
Xmax = max( x_element_index );
Ymin = min( y_element_index );
Ymax = max( y_element_index );
   
[ SAx, SAy ] = meshgrid( -Mx/2 : Mx/2 - 1, -My/2 : My/2 - 1 );

SA_array_x = ( Xmin + SAx + Mx/2 ) * spacing;
SA_array_y = ( Ymin + SAy + My/2 ) * spacing;

Rsam_SA = zeros( NelemSA );
cnt = 1;
   

recompute_OSA = 0;

if recompute_OSA
for mx = 1 : Nelem
    for nx = 1 : Nelem
        x_ind_SA = SAx + mx + Mx/2 - 1;
        y_ind_SA = SAy + nx + My/2 - 1;
        
        
        if( ( max( max( x_ind_SA + Xmin ) ) <= Xmax ) && ( min( min( x_ind_SA + Xmin ) ) >= Xmin ) && ( min( min( y_ind_SA + Ymin ) ) >= Ymin ) && ( max( max( y_ind_SA + Ymin ) ) <= Ymax ) )
        SA_index_X_store( cnt, 1 : Mx, 1 : My ) = x_ind_SA;
        SA_index_Y_store( cnt, 1 : Mx, 1 : My ) = y_ind_SA;

        SA_x_m( cnt, 1 : Mx, 1 : My ) = ( Xmin + x_ind_SA ) * spacing;
        SA_y_m( cnt, 1 : Mx, 1 : My ) = ( Ymin + y_ind_SA ) * spacing;

        figure( 555 )
        plot( array_x_planar, array_y_planar, 'bo' );
        hold on
        plot( ( Xmin + x_ind_SA ) * spacing, ( Ymin + y_ind_SA ) * spacing, 'ro' ); grid on
        hold off
        title('Overlapping Subarrays');
        xlabel('X'); ylabel('Y');
        pause( 0.01 )
        if ( print_plots && ( cnt == 75 ) )
            set( gcf, 'PaperPositionMode', 'auto' );
            eval( [ 'print( ''-dpng'', [ PlotPath, ''overlapped_subarrays.png''] )' ] );
            pause( 1 );
        end
             
        cnt = cnt + 1;
        end
    end
end
else
    display('Load the OSA layout');
    load osa_layout;
end

%save osa_layout SA_x_m SA_y_m SA_index_X_store SA_index_Y_store cnt


%%
NSAs = cnt - 1;
NSAx = Nelem / Mx;
NSAy = Nelem / My;
NSAs_non_overlap = NelemTot / ( Mx * My );


if steer_off == 1
    v_steer = v_steer_0 - 0.0152;
    u_steer = u_steer_0 - -0.0142;
elseif steer_off == 2
    v_steer = v_steer_0 - -0.0145;
    u_steer = u_steer_0 - 0.016;
elseif steer_off == 3
    v_steer = v_steer_0 - -0.0165;
    u_steer = u_steer_0 - -0.018;
elseif steer_off == 4
    v_steer = v_steer_0 - 0.0175;
    u_steer = u_steer_0 - 0.02;
else
    v_steer = v_steer_0;
    u_steer = u_steer_0;
end


display('Steer the array')
steering_phase = SA_steer( array_x_planar, array_y_planar, -u_steer, -v_steer, k );
steered_array_data0 = bsxfun( @times, steering_phase( : ).', Vt_new0 );

steered_array_data = Vt_new0;

Nfft = 256;
DeltaU_Nfft = [ -Nfft / 2 : Nfft / 2 - 1 ] * ( 2 / Nfft );  
DeltaV_Nfft = [ -Nfft / 2 : Nfft / 2 - 1 ] * ( 2 / Nfft );  

UplotFFT = zeros( Nfft, Nfft );
VplotFFT = zeros( Nfft, Nfft );
for lx = 1 : Nfft
    for nx = 1 : Nfft
        UplotFFT( lx, nx ) = DeltaU_Nfft( lx );
        VplotFFT( lx, nx ) = DeltaV_Nfft( nx );
    end
end
R2 = sqrt( UplotFFT .* UplotFFT + VplotFFT .* VplotFFT );

UplotFFT1 = UplotFFT;
VplotFFT1 = VplotFFT;

UplotFFT1( R2 > 1 ) = NaN;
VplotFFT1( R2 > 1 ) = NaN;


check_FTout = compute_FFT( UplotFFT1, VplotFFT1, squeeze( steered_array_data0( freq_sam, : ) ).', array_x_planar( : ), array_y_planar( : ), lambda );
check_FTout_dB = 20 * log10( abs( check_FTout ) );


%%
display('Pack steered array data into subarrays');
for ax = 1 : NSAs
    phase_center_X( ax, 1 ) = mean( mean( SA_x_m( ax, :, : ) ) );
    phase_center_Y( ax, 1 ) = mean( mean( SA_y_m( ax, :, : ) ) );
end


SA_data_steered = zeros( freq_sam, NSAs, Mx, My );
for tx = 1 : freq_sam
    Vt_temp = reshape( squeeze( steered_array_data( tx, 1 : NelemTot ) ), Nelem, Nelem );
    for ax = 1 : NSAs
        row_ind = squeeze( SA_index_X_store( ax, 1 : Mx, 1 : My ) );
        col_ind = squeeze( SA_index_Y_store( ax, 1 : Mx, 1 : My ) );
        for cx = 1 : Mx
            for dx = 1 : My
                SA_data_steered( tx, ax, cx, dx ) = Vt_temp( row_ind( cx, dx ) + 1, col_ind( cx, dx ) + 1 );
            end
        end
    end
end


%%
disp('Apply Errors')
Xmid = ( max( SA_array_x( : ) ) + min( SA_array_x( : ) ) ) / 2;
Ymid = ( max( SA_array_y( : ) ) + min( SA_array_y( : ) ) ) / 2;
distance_from_center_m = sqrt( ( SA_array_x - Xmid ).^2 + ( SA_array_y - Ymid ).^2 );

MaxPh_err_rad = pi/4;
MaxAmp_err_dB = 2;
MaxAmp_err = 10 ^ ( MaxAmp_err_dB / 20 );

Ph_err_rad = MaxPh_err_rad * sin( ( 2 * pi / lambda ) * ( distance_from_center_m / 10 ) );
Amp_err = MaxAmp_err * sin( ( 2 * pi / lambda ) * ( distance_from_center_m / 20 ) );

Ph_err_rad_2 = MaxPh_err_rad * rand( Mx, My );
Amp_err_2 = MaxAmp_err * randn( Mx, My );

if apply_errors == 1
    complexErrors = ( 1 + Amp_err ) .* exp( 1j *  Ph_err_rad );
elseif apply_errors == 2
    complexErrors = ( 1 + Amp_err_2 ) .* exp( 1j *  Ph_err_rad_2 );
else
    complexErrors = ones( Mx, My );
end

if applyBF_error == 1
    Ph_err_rad_3 = MaxPh_err_rad * rand;
    Amp_err_3 = MaxAmp_err * randn;
    complexBF_error = ( 1 + Amp_err_3 ) .* exp( 1j *  Ph_err_rad_3 );
else
    complexBF_error = 1;
end

SA_distance_from_center_m = sqrt( ( phase_center_X - XmidArray ).^2 + ( phase_center_Y - YmidArray ).^2 );

Ph_err_rad_SA = MaxPh_err_rad_3 * sin( ( 2 * pi / lambda ) * ( SA_distance_from_center_m / 15 ) );
Amp_err_SA = MaxAmp_err_3 * sin( ( 2 * pi / lambda ) * ( SA_distance_from_center_m / 30 ) );

if SA_phase_center_error == 1
    complexSAerrors = reshape( ( 1 + Amp_err_SA ) .* exp( 1j *  Ph_err_rad_SA ), NSAs, 1 );
else
    complexSAerrors = ones( NSAs, 1 );
end


SA_data_steered_out_vec = zeros( freq_sam, NSAs );
SA_steer_vec = SA_steer( phase_center_X, phase_center_Y, u_steer, v_steer, k );

for tx = 1 : freq_sam
    for ax = 1 : NSAs
        tempSA = squeeze( SA_data_steered( tx, ax, :, : ) );
        tempSA_errors = tempSA .* complexErrors;
        SA_data_steered_out_vec( tx, ax ) = SA_steer_vec( ax ) * sum( tempSA_errors( : ) ) * complexSAerrors( ax );
    end
end


check_FTout_SA = compute_FFT( UplotFFT1, VplotFFT1, squeeze( SA_data_steered_out_vec( freq_sam, : ) ).', phase_center_X, phase_center_Y, lambda );
check_FTout_SA_dB = 20 * log10( abs( check_FTout_SA ) );


if show_all_plots
   
figure( 1090 )
mesh( array_x_planar, array_y_planar, abs( reshape( complexPosErrors( 588, : ), Nelem, Nelem ) ) )
title('Position Errors')
xlabel('X [m]');
ylabel('Y [m]');


figure( 412 )
pcolor( UplotFFT1, VplotFFT1, check_FTout_dB );
xlabel( 'U' )
ylabel( 'V' );
title('Steered Array Output');
colormap jet; shading flat; axis('square'); colorbar; %caxis([ D11-65 D11 ])
if print_plots
    set( gcf, 'PaperPositionMode', 'auto' );
    eval( [ 'print( ''-dpng'', [ PlotPath, ''steered_array_output.png''] )' ] );
    pause( 1 );
end


figure( 492 )
pcolor( UplotFFT1, VplotFFT1, check_FTout_SA_dB );
xlabel( 'U' )
ylabel( 'V' );
title('Steered Subarray Output');
colormap jet; shading flat; axis('square'); colorbar; %caxis([ D11-65 D11 ])
if print_plots
    set( gcf, 'PaperPositionMode', 'auto' );
    eval( [ 'print( ''-dpng'', [ PlotPath, ''steered_subarray_output.png''] )' ] );
    pause( 1 );
end


figure( 493 )
plot( phase_center_X, phase_center_Y, 'o', 'MarkerSize', 4 );
title('Subarray Phase Centers');
xlabel( 'X [m]' );
ylabel( 'Y [m]' );
grid on
if print_plots
    set( gcf, 'PaperPositionMode', 'auto' );
    eval( [ 'print( ''-dpng'', [ PlotPath, ''subarray_phase_centers.png''] )' ] );
    pause( 1 );
end


figure( 494 )
mesh( reshape( phase_center_X, sqrt( NSAs ), sqrt( NSAs ) ), reshape( phase_center_Y, sqrt( NSAs ), sqrt( NSAs ) ), abs( reshape( complexSAerrors, sqrt( NSAs ), sqrt( NSAs ) ) ) );
title('Subarray Position Errors');
xlabel( 'X [m]' );
ylabel( 'Y [m]' );
grid on
if print_plots
    set( gcf, 'PaperPositionMode', 'auto' );
    eval( [ 'print( ''-dpng'', [ PlotPath, ''subarray_phase_center_errors.png''] )' ] );
    pause( 1 );
end

end         %  if show_all_plots


mainbeam_3dB = lambda / ( Nelem * lambda / 2 );      
mainbeam_3dB_0 = mainbeam_3dB; 

Nbeams = 16;  %8;

DeltaBeam = 2 * pi / Nbeams;
BeamSpacing = linspace( 0, 2 * pi - DeltaBeam, Nbeams );

beamU = cos( BeamSpacing );
beamV = sin( BeamSpacing );

beam_steer_ux = bsxfun( @plus, u_steer, mainbeam_3dB_0 * beamU.' );
beam_steer_vy = bsxfun( @plus, v_steer, mainbeam_3dB_0 * beamV.' );


%%  B MATRIX
display('Compute B matrix')

SA_vec = zeros( NSAs, Nbeams );
B = zeros( Nbeams, NSAs );

for ax = 1 : Nbeams
    tempS = ( 1 / sqrt( NSAs ) ) * SA_steer( beam_steer_ux( ax ), beam_steer_vy( ax ), phase_center_X, phase_center_Y, k );
    SA_vec( :, ax ) = tempS;
    %B( ax, : ) = tempS.' / ( tempS' * tempS );
    B( ax, : ) = tempS.';
end

%Rsam = SA_data_steered_out_vec.' * conj( SA_data_steered_out_vec );
Rsam = SA_data_steered_out_vec' * SA_data_steered_out_vec;

Snbr = 311;

SA_data_steered_after_IFT_0 = ifft( SA_data_steered_out_vec );
SA_data_steered_after_IFT_0_dB = 20 * log10( abs( SA_data_steered_after_IFT_0( :, Snbr ) ) );

SA_data_steered_after_IFT = SA_data_steered_after_IFT_0( :, Snbr );
[ I1, I2 ] = max( 20 * log10( abs( SA_data_steered_after_IFT ) ) );

Bin_nbr = I2
xi = B * SA_data_steered_after_IFT_0.';
psi_0 = xi( :, Bin_nbr );
psi_0_dB = 20 * log10( abs( psi_0 ) );


if show_all_plots
    
figure( 458 )
plot( psi_0_dB, 'LineWidth', 2 );
grid on
xlabel( 'Beam Index' )
ylabel( 'dB' );
axis([ 1 Nbeams min( psi_0_dB )-3 max( psi_0_dB )+3 ])
title('Beam Outputs');
if print_plots
    set( gcf, 'PaperPositionMode', 'auto' );
    eval( [ 'print( ''-dpng'', [ PlotPath, ''psi_0.png''] )' ] );
    pause( 1 );
end


figure( 460 )
plot( SA_data_steered_after_IFT_0_dB, 'LineWidth', 1.6 );
grid on
xlabel( 'Delay Index' )
ylabel( 'dB' );
title('Subarray Output');
axis([ -1 freq_sam min( SA_data_steered_after_IFT_0_dB )-5 max( SA_data_steered_after_IFT_0_dB )+5 ])
if print_plots
    set( gcf, 'PaperPositionMode', 'auto' );
    eval( [ 'print( ''-dpng'', [ PlotPath, ''SA_ifft_output.png''] )' ] );
    pause( 1 );
end


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

end



%%
display('Start Beam of Beams');

C = B * Rsam * B';
%C = B * B';

rank_C = rank( C )
rank_B = rank( B )
rank_Rsam = rank( Rsam )

nu_3 = 0;
C_inv = inv( C + nu_3 * eye( Nbeams ) );


%%
Duv = 0.5;
Umin = -Duv;  
Vmin = -Duv; 
Umax = Duv;  
Vmax = Duv;  

Ndelta = 100;
cost_func = zeros( Ndelta, Ndelta );
ugrid = linspace( Umin, Umax, Ndelta );
vgrid = linspace( Vmin, Vmax, Ndelta );

[ DeltaUgrid, DeltaVgrid ] = meshgrid(  ugrid, vgrid  );


%% cost_function
display('Compute cost function');

init_phases = zeros( 1, NSAs );

for jx = 1 : Ndelta
    for lx = 1 : Ndelta
        uv0 = [ u_steer + DeltaUgrid( jx, lx ), v_steer + DeltaVgrid( jx, lx ) ];
        cost_func( jx, lx ) = SA_cost_function_SteerVec_plot( uv0, k, B, C_inv, psi_0, u_steer, v_steer, phase_center_X, phase_center_Y );
    end
end
cost_func_dB = 20 * log10( abs( cost_func ) );



%%
display('Create zoom plot');
Duv2 = 0.04;  
Umin2 = -Duv2;  
Vmin2 = -Duv2; 
Umax2 = Duv2;  
Vmax2 = Duv2;  

Ndelta2 = 100;
cost_func2 = zeros( Ndelta2, Ndelta2 );
ugrid2 = linspace( Umin2, Umax2, Ndelta2 );
vgrid2 = linspace( Vmin2, Vmax2, Ndelta2 );

[ DeltaUgrid2, DeltaVgrid2 ] = meshgrid(  ugrid2, vgrid2  );


%% cost_function
display('Compute cost function');

for jx = 1 : Ndelta2
    for lx = 1 : Ndelta2
        uv00 = [ u_steer + DeltaUgrid2( jx, lx ), v_steer + DeltaVgrid2( jx, lx ) ];
        cost_func2( jx, lx ) = SA_cost_function_SteerVec_plot( uv00, k, B, C_inv, psi_0, u_steer, v_steer, phase_center_X, phase_center_Y );
    end
end
cost_func_dB2 = 20 * log10( abs( cost_func2 ) );


%% search for the optimum
options = optimset( 'MaxIter', 2e2, 'Display', 'iter', 'GradObj', 'off', 'Hessian', 'off', 'DerivativeCheck', 'off', 'TolFun', 1e-16, 'TolX', 1e-16, 'MaxFunEvals', 2e5 );
options2 = optimset( 'MaxIter', 200, 'Display', 'iter', 'GradObj', 'off', 'Hessian', 'off', 'DerivativeCheck', 'off', 'TolFun', 1e-40, 'TolX', 1e-40, 'MaxFunEvals', 2e5 );

SA_cost_function_1 = @( x )SA_cost_function_SteerVec( x, k, B, C_inv, psi_0, u_steer, v_steer, phase_center_X, phase_center_Y );

[ PkCost, Idx ] = max( cost_func_dB2( : ) );
[ r0, c0 ] = ind2sub( [ Ndelta2, Ndelta2 ], Idx );
u00 = u_steer + DeltaUgrid2( r0, c0 );
v00 = v_steer + DeltaVgrid2( r0, c0 );

x00 = [ u_steer, v_steer ];
%x00 = [ u00, v00 ];

[ ML, psi_out, exitflag ] = fminsearch( SA_cost_function_1, x00, options );

u_est = ML( 1 ) 
v_est = ML( 2 ) 


%%
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




if show_all_plots
%%
figure( 122 )
plot3( u_steer + DeltaUgrid, v_steer + DeltaVgrid, cost_func_dB )
grid on
xlabel( 'U' ); ylabel( 'V' ); zlabel( 'Cost (dB)' );
title( 'Cost' )
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
title( 'PColor Plot of Cost' )
colormap jet; colorbar; caxis([ L1 L2 ]);  %caxis([ min( cost_func_dB( : ) ) max( cost_func_dB( : ) ) ])
if print_plots
    set( gcf, 'PaperPositionMode', 'auto' );
    eval( [ 'print( ''-dpng'', [ PlotPath, ''cost_mle_pcolor.png''] )' ] );
    pause( 1 )
end


%%
L11 = min( cost_func_dB2( : ) ) - 5;
L22 = max( cost_func_dB2( : ) ) + 5;
figure( 125 )
pcolor( u_steer + DeltaUgrid2, v_steer + DeltaVgrid2, cost_func_dB2 )
grid on
xlabel( 'U' ); ylabel( 'V' ); zlabel( 'Cost (dB)' );
title( 'Zoomed Plot in Pcolor' )
colormap jet; colorbar; caxis([ L11 L22 ]);  %caxis([ min( cost_func_dB( : ) ) max( cost_func_dB( : ) ) ])
% view(90,0)
if print_plots
    set( gcf, 'PaperPositionMode', 'auto' );
    eval( [ 'print( ''-dpng'', [ PlotPath, ''cost_mle_pcolor_v2.png''] )' ] );
    pause( 1 )
end


%%
figure( 124 )
plot3( u_steer + DeltaUgrid2, v_steer + DeltaVgrid2, cost_func_dB2 )
grid on
xlabel( 'U' ); ylabel( 'V' ); zlabel( 'Cost (dB)' );
title( '3D Zoomed Plot of Cost' )
view( 62, 80 )
if print_plots
    set( gcf, 'PaperPositionMode', 'auto' );
    eval( [ 'print( ''-dpng'', [ PlotPath, ''cost_mle_v2.png''] )' ] );
    pause( 1 )
end


end


figure( 126 )
plot( beam_steer_ux, beam_steer_vy, 'o', 'MarkerSize', 7, 'MarkerFace', [ 0.02 0.36 0.49 ] );
hold on
plot( u_steer, v_steer, 'rx', 'MarkerSize', 7 )
plot( u_steer_0, v_steer_0, 'kx', 'MarkerSize', 7 );
plot( u_est, v_est, 'bx', 'MarkerSize', 7 );
hold off
grid on
xlabel('U');
ylabel('V');
eval( [ 'title( ''Rosette Beam Locations:  Bandwidth = ', num2str( BandWidth_Hz / 1e9 ), ' GHz'', ''FontSize'', 14, ''Interpreter'', ''Latex'' )' ] );
legend('Rosette', 'Steer', 'Object', 'Estimate' );
if print_plots
    set( gcf, 'PaperPositionMode', 'auto' );
    eval( [ 'print( ''-dpng'', [ PlotPath, ''rosette_beams_', num2str( BandWidth_Hz / 1e9 ), '.png''] )' ] );
    pause( 3 )
end


if show_all_plots_2
figure( 1241 )
plot3( u_steer + DeltaUgrid2, v_steer + DeltaVgrid2, cost_func_dB2 )
grid on
xlabel( 'U' ); ylabel( 'V' ); zlabel( 'Cost (dB)' );
title( 'Zoomed Plot of Cost along the U-Axis' )
view(90,0)
if print_plots
    set( gcf, 'PaperPositionMode', 'auto' );
    eval( [ 'print( ''-dpng'', [ PlotPath, ''cost_mle_v2_view_90_0.png''] )' ] );
    pause( 1 )
end


figure( 1242 )
plot3( u_steer + DeltaUgrid2, v_steer + DeltaVgrid2, cost_func_dB2 )
grid on
xlabel( 'U' ); ylabel( 'V' ); zlabel( 'Cost (dB)' );
title( 'Zoomed Plot of Cost along the V-Axis' )
view(0,0)
if print_plots
    set( gcf, 'PaperPositionMode', 'auto' );
    eval( [ 'print( ''-dpng'', [ PlotPath, ''cost_mle_v2_view_0_90.png''] )' ] );
    pause( 1 )
end
end



%save CASE_4_OUTPUT_041620







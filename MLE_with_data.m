%%  This script feeds measured data through the MLE beamspace algorithm
clc
clear all
close all


print_plots = 0;

f_Hz = 40e9;  %28e9;  %1e9;
c = 3e8;
lambda = c / f_Hz;
k = 2 * pi / lambda;
spacing = lambda/2;


[ Vt_new, array_x_planar, array_y_planar, array_z_planar ] = save_Samurai_data_file;


freq_sam = 1350;
Nelem = 35;
NelemTot = Nelem * Nelem;


Nx = 7;
Ny = 7;
NelemSA = Nx * Ny;

NSAx = Nelem / Nx;
NSAy = Nelem / Ny;
NSAs = NelemTot / ( Nx * Ny );

SA_steer = @(x,y,u,v,k)( exp(-1j*k*(x.*u + y.*v)) );

u_steer = 0.0571;       %  these values correspond to peak
v_steer = 0.4857;

display('Steer the array')
steering_phase = SA_steer( array_x_planar, array_y_planar, u_steer, v_steer, k );
steered_array_data = bsxfun( @times, steering_phase( : )', Vt_new );

SA_ind_R = zeros( NSAs, Nx );
SA_ind_C = zeros( NSAs, Ny );
SA_data_steered = zeros( NSAs, Nx, Ny, freq_sam );

for tx = 1 : freq_sam
    Vt_temp = reshape( squeeze( steered_array_data( tx, 1 : NelemTot ) ), Nelem, Nelem );
    for ax = 1 : NSAx
        for bx = 1 : NSAy
        row_ind = ( ax - 1 ) * Nx + 1 : ax * Nx;
        col_ind = ( bx - 1 ) * Ny + 1 : bx * Ny;
        SA_ind_R( ( ax - 1 ) * NSAy + bx, : ) = row_ind;
        SA_ind_C( ( ax - 1 ) * NSAy + bx, : ) = col_ind;
        SA_x_m( ( ax - 1 ) * NSAy + bx, : ) = reshape( array_x_planar( row_ind, col_ind ), 1, NelemSA );
        SA_y_m( ( ax - 1 ) * NSAy + bx, : ) = reshape( array_y_planar( row_ind, col_ind ), 1, NelemSA );
        SA_data_steered( ( ax - 1 ) * NSAy + bx, 1 : Nx, 1 : Ny, tx ) = Vt_temp( row_ind, col_ind );
        end
    end
end


%%
show_array = 1;
if show_array
figure(998)
plot( array_x_planar, array_y_planar, 'ok' ); grid on
end


SA_data_steered_out_vec = zeros( freq_sam, NSAs );
for tx = 1 : freq_sam
    for ax = 1 : NSAs
        tempSA = squeeze( SA_data_steered( ax, :, :, tx ) );
        SA_data_steered_out_vec( tx, ax ) = sum( tempSA( : ) );
    end
end


mainbeam_3dB = lambda / ( Nelem * lambda );      
mainbeam_3dB_0 = mainbeam_3dB / 2;  
beam_steer_ux = bsxfun( @plus, u_steer, [ mainbeam_3dB_0 0 -mainbeam_3dB_0 0 ].');
beam_steer_vy = bsxfun( @plus, v_steer, [ 0 mainbeam_3dB_0 0 -mainbeam_3dB_0 ].');
 
%%  B MATRIX
display('Compute B matrix')
Nbeams = 4;

SA_vec = zeros( NSAs, Nbeams );
B = zeros( Nbeams, NSAs );
for ax = 1 : Nbeams
    for bx = 1 : NSAs
        tempX = squeeze( SA_x_m( bx, : ) );
        tempY = squeeze( SA_y_m( bx, : ) );
        SA_vec( bx, ax ) = ( 1 / sqrt( NelemSA ) ) * sum( SA_steer( tempX, tempY, beam_steer_ux( ax ) - u_steer, beam_steer_vy( ax ) - v_steer, k ) );
    end
end
        
for ax = 1 : Nbeams
    tempS = squeeze( SA_vec( 1 : NSAs, ax ) );
    B( ax, : ) = tempS / ( tempS' * tempS );
end

Rsam = SA_data_steered_out_vec.' * conj( SA_data_steered_out_vec );

Snbr = 311;
xi = B * SA_data_steered_out_vec.';
psi_0 = xi( :, Snbr );

C = B * Rsam * B';

rank_C = rank( C )
rank_B = rank( B )

%C_inv = eye( Nbeams );
nu_3 = 0;
C_inv = inv( C + nu_3 * eye( Nbeams ) );

%% loop in Duv

% Duv_loop = [0.01:0.01:0.5];  
% for i = 1: length(Duv_loop)
% Duv =  Duv_loop(i)   ; %0.04;  %0.5;
% Umin = -Duv;  
% Vmin = -Duv; 
% Umax = Duv;  
% Vmax = Duv;  
% 
% Ndelta = 200; 
% cost_func = zeros( Ndelta, Ndelta );
% ugrid = linspace( Umin, Umax, Ndelta );
% vgrid = linspace( Vmin, Vmax, Ndelta );
% 
% [ DeltaUgrid, DeltaVgrid ] = meshgrid(  ugrid, vgrid  );
% 
% % end
% %% cost_function onuvgrid 
% display('Compute cost function');
% 
% for jx = 1 : Ndelta
%     for kx = 1 : Ndelta
%         uv0 = [ u_steer + DeltaUgrid( jx, kx ), v_steer + DeltaVgrid( jx, kx ) ];
%         cost_func( jx, kx ) = 
%SA_cost_function_data( uv0, k, NSAs, B, C_inv, SA_x_m, SA_y_m, psi_0, u_steer, v_steer );
%     end
% end
% cost_func_dB = 20 * log10( abs( cost_func ) );
% 
% PlotPath = '\\bd5\5g-share\5G_Share\Abhishek\E_band_feb12\temp\MVDR\Matlab_Code_for_Abhishek_091219\fig_beamspace_real_data\Samurai_Conference_Room\' ;
% print_plots = 1;
% 
% figure( 122 )
% plot3( u_steer + DeltaUgrid, v_steer + DeltaVgrid, cost_func_dB )
% grid on
% xlabel( 'U', 'FontSize', 14, 'Interpreter', 'Latex' ); ylabel( 'V', 'FontSize', 14, 'Interpreter', 'Latex' ); zlabel( 'Cost (dB)' , 'FontSize', 14, 'Interpreter', 'Latex' );
% eval( [ 'title( ''Cost Function in Steering Grid Using plot3 Duv = ', num2str(Duv) , ''', ''FontSize'', 14, ''Interpreter'', ''Latex'' )' ] );
% view( 34, 17 )
% 
% if print_plots
%     set( gcf, 'PaperPositionMode', 'auto' );
%     eval( [ 'print( ''-dpng'', [ PlotPath, ''cost_varying_grid_', num2str( i ), '.png''] )' ] ); pause( 1 )
% end
% 
% %%
% figure( 123 )
% pcolor( u_steer + DeltaUgrid, v_steer + DeltaVgrid, cost_func_dB )
% grid on
% xlabel( 'U', 'FontSize', 14, 'Interpreter', 'Latex' ); ylabel( 'V', 'FontSize', 14, 'Interpreter', 'Latex' ); zlabel( 'Cost (dB)', 'FontSize', 14, 'Interpreter', 'Latex' );
% eval( [ 'title( ''Cost Function in Steering Grid Using pcolor Duv = ', num2str( Duv ) , ''',''FontSize'', 14, ''Interpreter'', ''Latex'' )' ] );
% colormap jet; colorbar; caxis([ 70 150 ]);  %caxis([ min( cost_func_dB( : ) ) max( cost_func_dB( : ) ) ])
% 
% if print_plots
%     set( gcf, 'PaperPositionMode', 'auto' );
%     eval( [ 'print( ''-dpng'', [ PlotPath, ''cost_varying_grid_pcolor_', num2str( i ), '.png''] )' ] ); pause( 1 )
% end
% 
% end
%% obs
fig_vallid_Duv = [0.01:0.01:0.05];  %opt performance range for conference room data bob real 

Duv_loop = [0.04];          %best value of opt range   

u_steer_loop = [ 0.051 : 0.001 : 0.060];
v_steer 

for j = 1: length(u_steer_loop)
Duv =  Duv_loop   ; %0.04;  %0.5;
Umin = -Duv;  
Vmin = -Duv; 
Umax = Duv;  
Vmax = Duv;  

Ndelta = 200; 
cost_func = zeros( Ndelta, Ndelta );
ugrid = linspace( Umin, Umax, Ndelta );
vgrid = linspace( Vmin, Vmax, Ndelta );

[ DeltaUgrid, DeltaVgrid ] = meshgrid(  ugrid, vgrid  );

u_steer = u_steer_loop(j) ;

%  end
%% cost_function onuvgrid 
display('Compute cost function');

for jx = 1 : Ndelta
    for kx = 1 : Ndelta
        uv0 = [ u_steer + DeltaUgrid( jx, kx ), v_steer + DeltaVgrid( jx, kx ) ];
        cost_func( jx, kx ) = SA_cost_function_data( uv0, k, NSAs, B, C_inv, SA_x_m, SA_y_m, psi_0, u_steer, v_steer );
    end
end
cost_func_dB = 20 * log10( abs( cost_func ) );

PlotPath = '\\bd5\5g-share\5G_Share\Abhishek\E_band_feb12\temp\MVDR\Matlab_Code_for_Abhishek_091219\fig_beamspace_real_data\Samurai_Conference_Room\' ;
print_plots = 1;

figure( 122 )
plot3( u_steer + DeltaUgrid, v_steer + DeltaVgrid, cost_func_dB )
grid on
xlabel( 'U', 'FontSize', 14, 'Interpreter', 'Latex' ); ylabel( 'V', 'FontSize', 14, 'Interpreter', 'Latex' ); zlabel( 'Cost (dB)' , 'FontSize', 14, 'Interpreter', 'Latex' );
eval( [ 'title( ''Cost Function in Steering Grid Using plot3 u_steer = ', num2str(u_steer) , ''', ''FontSize'', 14, ''Interpreter'', ''Latex'' )' ] );
view( 34, 17 )

if print_plots
    set( gcf, 'PaperPositionMode', 'auto' );
    eval( [ 'print( ''-dpng'', [ PlotPath, ''cost_varying_u_steer_', num2str( j ), '.png''] )' ] ); pause( 1 )
end

%%
figure( 123 )
pcolor( u_steer + DeltaUgrid, v_steer + DeltaVgrid, cost_func_dB )
grid on
xlabel( 'U', 'FontSize', 14, 'Interpreter', 'Latex' ); ylabel( 'V', 'FontSize', 14, 'Interpreter', 'Latex' ); zlabel( 'Cost (dB)', 'FontSize', 14, 'Interpreter', 'Latex' );
eval( [ 'title( ''Cost Function in Steering Grid Using pcolor u_steer = ', num2str( u_steer ) , ''',''FontSize'', 14, ''Interpreter'', ''Latex'' )' ] );
colormap jet; colorbar; caxis([ 70 150 ]);  %caxis([ min( cost_func_dB( : ) ) max( cost_func_dB( : ) ) ])

if print_plots
    set( gcf, 'PaperPositionMode', 'auto' );
    eval( [ 'print( ''-dpng'', [ PlotPath, ''cost_varying_u_steer_pcolor_', num2str( j ), '.png''] )' ] ); pause( 1 )
end

end













print_plots = 0;

f_Hz = 40e9;  %28e9;  %1e9;
c = 3e8;
lambda = c / f_Hz;
k = 2 * pi / lambda;
spacing = lambda/2;


[ Vt_new, array_x_planar, array_y_planar, array_z_planar ] = save_Samurai_data_file;


freq_sam = 1350;
Nelem = 35;
NelemTot = Nelem * Nelem;


Nx = 7;
Ny = 7;
NelemSA = Nx * Ny;

NSAx = Nelem / Nx;
NSAy = Nelem / Ny;
NSAs = NelemTot / ( Nx * Ny );

SA_steer = @(x,y,u,v,k)( exp(-1j*k*(x.*u + y.*v)) );

u_steer = 0.0571;       %  these values correspond to peak
v_steer = 0.4857;

display('Steer the array')
steering_phase = SA_steer( array_x_planar, array_y_planar, u_steer, v_steer, k );
steered_array_data = bsxfun( @times, steering_phase( : )', Vt_new );

SA_ind_R = zeros( NSAs, Nx );
SA_ind_C = zeros( NSAs, Ny );
SA_data_steered = zeros( NSAs, Nx, Ny, freq_sam );

for tx = 1 : freq_sam
    Vt_temp = reshape( squeeze( steered_array_data( tx, 1 : NelemTot ) ), Nelem, Nelem );
    for ax = 1 : NSAx
        for bx = 1 : NSAy
        row_ind = ( ax - 1 ) * Nx + 1 : ax * Nx;
        col_ind = ( bx - 1 ) * Ny + 1 : bx * Ny;
        SA_ind_R( ( ax - 1 ) * NSAy + bx, : ) = row_ind;
        SA_ind_C( ( ax - 1 ) * NSAy + bx, : ) = col_ind;
        SA_x_m( ( ax - 1 ) * NSAy + bx, : ) = reshape( array_x_planar( row_ind, col_ind ), 1, NelemSA );
        SA_y_m( ( ax - 1 ) * NSAy + bx, : ) = reshape( array_y_planar( row_ind, col_ind ), 1, NelemSA );
        SA_data_steered( ( ax - 1 ) * NSAy + bx, 1 : Nx, 1 : Ny, tx ) = Vt_temp( row_ind, col_ind );
        end
    end
end


%%
show_array = 1;
if show_array
figure(998)
plot( array_x_planar, array_y_planar, 'ok' ); grid on
end


SA_data_steered_out_vec = zeros( freq_sam, NSAs );
for tx = 1 : freq_sam
    for ax = 1 : NSAs
        tempSA = squeeze( SA_data_steered( ax, :, :, tx ) );
        SA_data_steered_out_vec( tx, ax ) = sum( tempSA( : ) );
    end
end


mainbeam_3dB = lambda / ( Nelem * lambda );      
mainbeam_3dB_0 = mainbeam_3dB / 2;  
beam_steer_ux = bsxfun( @plus, u_steer, [ mainbeam_3dB_0 0 -mainbeam_3dB_0 0 ].');
beam_steer_vy = bsxfun( @plus, v_steer, [ 0 mainbeam_3dB_0 0 -mainbeam_3dB_0 ].');
 


%%  B MATRIX
display('Compute B matrix')
Nbeams = 4;

SA_vec = zeros( NSAs, Nbeams );
B = zeros( Nbeams, NSAs );
for ax = 1 : Nbeams
    for bx = 1 : NSAs
        tempX = squeeze( SA_x_m( bx, : ) );
        tempY = squeeze( SA_y_m( bx, : ) );
        SA_vec( bx, ax ) = ( 1 / sqrt( NelemSA ) ) * sum( SA_steer( tempX, tempY, beam_steer_ux( ax ) - u_steer, beam_steer_vy( ax ) - v_steer, k ) );
    end
end
        
for ax = 1 : Nbeams
    tempS = squeeze( SA_vec( 1 : NSAs, ax ) );
    B( ax, : ) = tempS / ( tempS' * tempS );
end

Rsam = SA_data_steered_out_vec.' * conj( SA_data_steered_out_vec );

Snbr = 311;
xi = B * SA_data_steered_out_vec.';
psi_0 = xi( :, Snbr );

C = B * Rsam * B';

rank_C = rank( C )
rank_B = rank( B )

%C_inv = eye( Nbeams );
nu_3 = 0;
C_inv = inv( C + nu_3 * eye( Nbeams ) );

%% loop in Duv

% Duv_loop = [0.01:0.01:0.5];  
% for i = 1: length(Duv_loop)
% Duv =  Duv_loop(i)   ; %0.04;  %0.5;
% Umin = -Duv;  
% Vmin = -Duv; 
% Umax = Duv;  
% Vmax = Duv;  
% 
% Ndelta = 200; 
% cost_func = zeros( Ndelta, Ndelta );
% ugrid = linspace( Umin, Umax, Ndelta );
% vgrid = linspace( Vmin, Vmax, Ndelta );
% 
% [ DeltaUgrid, DeltaVgrid ] = meshgrid(  ugrid, vgrid  );
% 
% % end
% %% cost_function onuvgrid 
% display('Compute cost function');
% 
% for jx = 1 : Ndelta
%     for kx = 1 : Ndelta
%         uv0 = [ u_steer + DeltaUgrid( jx, kx ), v_steer + DeltaVgrid( jx, kx ) ];
%         cost_func( jx, kx ) = 
%SA_cost_function_data( uv0, k, NSAs, B, C_inv, SA_x_m, SA_y_m, psi_0, u_steer, v_steer );
%     end
% end
% cost_func_dB = 20 * log10( abs( cost_func ) );
% 
% PlotPath = '\\bd5\5g-share\5G_Share\Abhishek\E_band_feb12\temp\MVDR\Matlab_Code_for_Abhishek_091219\fig_beamspace_real_data\Samurai_Conference_Room\' ;
% print_plots = 1;
% 
% figure( 122 )
% plot3( u_steer + DeltaUgrid, v_steer + DeltaVgrid, cost_func_dB )
% grid on
% xlabel( 'U', 'FontSize', 14, 'Interpreter', 'Latex' ); ylabel( 'V', 'FontSize', 14, 'Interpreter', 'Latex' ); zlabel( 'Cost (dB)' , 'FontSize', 14, 'Interpreter', 'Latex' );
% eval( [ 'title( ''Cost Function in Steering Grid Using plot3 Duv = ', num2str(Duv) , ''', ''FontSize'', 14, ''Interpreter'', ''Latex'' )' ] );
% view( 34, 17 )
% 
% if print_plots
%     set( gcf, 'PaperPositionMode', 'auto' );
%     eval( [ 'print( ''-dpng'', [ PlotPath, ''cost_varying_grid_', num2str( i ), '.png''] )' ] ); pause( 1 )
% end
% 
% %%
% figure( 123 )
% pcolor( u_steer + DeltaUgrid, v_steer + DeltaVgrid, cost_func_dB )
% grid on
% xlabel( 'U', 'FontSize', 14, 'Interpreter', 'Latex' ); ylabel( 'V', 'FontSize', 14, 'Interpreter', 'Latex' ); zlabel( 'Cost (dB)', 'FontSize', 14, 'Interpreter', 'Latex' );
% eval( [ 'title( ''Cost Function in Steering Grid Using pcolor Duv = ', num2str( Duv ) , ''',''FontSize'', 14, ''Interpreter'', ''Latex'' )' ] );
% colormap jet; colorbar; caxis([ 70 150 ]);  %caxis([ min( cost_func_dB( : ) ) max( cost_func_dB( : ) ) ])
% 
% if print_plots
%     set( gcf, 'PaperPositionMode', 'auto' );
%     eval( [ 'print( ''-dpng'', [ PlotPath, ''cost_varying_grid_pcolor_', num2str( i ), '.png''] )' ] ); pause( 1 )
% end
% 
% end
%% obs
fig_vallid_Duv = [0.01:0.01:0.05];  %opt performance range for conference room data bob real 

Duv_loop = [0.04];          %best value of opt range   

u_steer_loop = [ 0.051 : 0.001 : 0.060];
v_steer 

for j = 1: length(u_steer_loop)
Duv =  Duv_loop   ; %0.04;  %0.5;
Umin = -Duv;  
Vmin = -Duv; 
Umax = Duv;  
Vmax = Duv;  

Ndelta = 200; 
cost_func = zeros( Ndelta, Ndelta );
ugrid = linspace( Umin, Umax, Ndelta );
vgrid = linspace( Vmin, Vmax, Ndelta );

[ DeltaUgrid, DeltaVgrid ] = meshgrid(  ugrid, vgrid  );

u_steer = u_steer_loop(j) ;

%  end
%% cost_function onuvgrid 
display('Compute cost function');

for jx = 1 : Ndelta
    for kx = 1 : Ndelta
        uv0 = [ u_steer + DeltaUgrid( jx, kx ), v_steer + DeltaVgrid( jx, kx ) ];
        cost_func( jx, kx ) = SA_cost_function_data( uv0, k, NSAs, B, C_inv, SA_x_m, SA_y_m, psi_0, u_steer, v_steer );
    end
end
cost_func_dB = 20 * log10( abs( cost_func ) );

PlotPath = '\\bd5\5g-share\5G_Share\Abhishek\E_band_feb12\temp\MVDR\Matlab_Code_for_Abhishek_091219\fig_beamspace_real_data\Samurai_Conference_Room\' ;
print_plots = 1;

figure( 122 )
plot3( u_steer + DeltaUgrid, v_steer + DeltaVgrid, cost_func_dB )
grid on
xlabel( 'U', 'FontSize', 14, 'Interpreter', 'Latex' ); ylabel( 'V', 'FontSize', 14, 'Interpreter', 'Latex' ); zlabel( 'Cost (dB)' , 'FontSize', 14, 'Interpreter', 'Latex' );
eval( [ 'title( ''Cost Function in Steering Grid Using plot3 u_steer = ', num2str(u_steer) , ''', ''FontSize'', 14, ''Interpreter'', ''Latex'' )' ] );
view( 34, 17 )

if print_plots
    set( gcf, 'PaperPositionMode', 'auto' );
    eval( [ 'print( ''-dpng'', [ PlotPath, ''cost_varying_u_steer_', num2str( j ), '.png''] )' ] ); pause( 1 )
end

%%
figure( 123 )
pcolor( u_steer + DeltaUgrid, v_steer + DeltaVgrid, cost_func_dB )
grid on
xlabel( 'U', 'FontSize', 14, 'Interpreter', 'Latex' ); ylabel( 'V', 'FontSize', 14, 'Interpreter', 'Latex' ); zlabel( 'Cost (dB)', 'FontSize', 14, 'Interpreter', 'Latex' );
eval( [ 'title( ''Cost Function in Steering Grid Using pcolor u_steer = ', num2str( u_steer ) , ''',''FontSize'', 14, ''Interpreter'', ''Latex'' )' ] );
colormap jet; colorbar; caxis([ 70 150 ]);  %caxis([ min( cost_func_dB( : ) ) max( cost_func_dB( : ) ) ])

if print_plots
    set( gcf, 'PaperPositionMode', 'auto' );
    eval( [ 'print( ''-dpng'', [ PlotPath, ''cost_varying_u_steer_pcolor_', num2str( j ), '.png''] )' ] ); pause( 1 )
end

end












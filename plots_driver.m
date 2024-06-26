%%  This script creates a plot of AOA error vs signal bandwidth




clear all
close all
clc


%unix('cp -Rf *.m BACKUP');

print_plots = 0;

PlotPath = '/home/pgv/Davis_Fante_Simulation/With_Position_Errors/Fast_Version/Plots_vs_Bandwidth_1_GHz_CASE_3_041620/';

%BandWidth_Hz = 0 : 1e9 : 14e9;
%Nf = length( BandWidth_Hz );

Nf = 15;
BandWidth_Hz = linspace( 0, 1e9, Nf );


for ax = 1 : Nf
    [ rms_error( ax ), rms_error_0( ax ) ] = fast_beamspace_MLE_version( BandWidth_Hz( ax ), PlotPath );
end



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








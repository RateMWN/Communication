function FTout = compute_FFT( U, V, temp, array_x, array_y, lambda )




F1 = U / lambda;
F2 = V / lambda;

[ NelemX, NelemY ] = size( array_x );
[ Nu, Nv ] = size( U );

Npos = NelemX * NelemY;


%%  Compute the FFT2
FTout = zeros( Nu, Nv );

for jx = 1 : NelemX
    for kx = 1 : NelemY
        temp2D = temp( jx, kx ) * exp( -1j * 2 * pi * ( array_x( jx, kx ) * F1 + array_y( jx, kx ) * F2 ) );
        FTout = FTout + temp2D;
    end
end


FTout = ( 1 / sqrt( Npos ) ) * FTout;



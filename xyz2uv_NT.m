function [ u, v ] = xyz2uv_NT( x, y, z, N, T )


R = sqrt( x * x + y * y + z * z );

u = ( x / R ) * cos( N ) - ( z / R ) * sin( N );

v = -( x / R ) * sin( N ) * sin( T ) + ( y / R ) * cos( T ) - ( z / R ) * cos( N ) * sin( T );

end

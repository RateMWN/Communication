function [ x, y, z ] = Razel2xyz( R, az, el )

x = R .* cos( el ) .* sin( az );
y = R .* sin( el );
z = R .* cos( el ) .* cos( az );

end

function [ u, v ] = azel2uv( az, el )

u = cos( el ) .* sin( az );
v = sin( el );

end

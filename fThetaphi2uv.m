function [u v] = fThetaphi2uv(theta, phi)
  %%
   N = length(phi);
   u = sin(theta).*cos(phi);
   v = sin(theta).*sin(phi);
   z_uv = cos(theta);%No need of z in uv because these are in  array surface
%    fprintf('fThetaphi2uv op is %d.\n',N);
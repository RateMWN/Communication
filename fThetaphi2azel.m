%version change
   function [az el] = fThetaphi2azel(theta,phi)
       az = atan(tan(theta).*cos(phi));    
       el = asin(sin(theta).*sin(phi)) ;
    function [x y z] = fThetaphi2xy(theta, phi)
      global r
           x = r.*sin(theta).*cos(phi);
           y = r.*sin(theta).*sin(phi);
           z = sqrt(r.^2-(x.^2+y.^2));
    

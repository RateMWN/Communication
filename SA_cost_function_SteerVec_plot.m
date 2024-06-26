function cost = SA_cost_function_SteerVec_plot( thPhi, k_wn, Bn, Cn_inv, psi_0, u_steer, v_steer, SA_phase_center_X, SA_phase_center_Y )


SAphase = @(u,v,x,y,k)exp(-1j*k*(x*u+y*v));


search_ux = -u_steer + thPhi( 1 );              
search_vy = -v_steer + thPhi( 2 );

steer_vec = SAphase( search_ux, search_vy, SA_phase_center_X, SA_phase_center_Y, k_wn );

Gn = Bn * steer_vec;

denom = ( Gn' * Cn_inv * Gn );

Tn = ( Cn_inv * Gn ) / denom;

SIR = denom * abs( Tn' * psi_0 ) ^ 2;

%cost = 1 / real( SIR );
cost = real( SIR );

end



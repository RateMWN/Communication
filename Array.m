%% plot all the SA in x-y plane
plot_array = 1;
%plot_array = 0;
%if plot_array
for short = 1:1
    
origin_vec =  [0: Nx-1]
origin_pos =  [-1 0 0; 0 0 0; 0 0 0]
origin_next_pos = eye (NSAs_x) +  origin_pos;
% 
initial_index = Nx * origin_next_pos * [0: NSAs_x-1 ].';
% next_increment= repmat( initial_index, [1, Nx]  );
% 
% order_xy = repmat(origin_vec, [NSAs_x, 1])+ next_increment;
% order_SA = order_xy* spacing;
% 
% [xall_SA yall_SA] = meshgrid( order_SA, order_SA);
end %replacementof last end
% %% plot max min cent0 of SA
% 
% atemp = [1:Nx];         %initial index, it starts from 1 not 0
% ia = NSAs ;
% SA_x0 =zeros (NSAs,1) ;  
% SA_y0 = zeros(NSAs,1);  array_x0= zeros(NSAs) ; array_y0 =zeros (NSAs);
% %%
%  for ai= 1: NSAs         % vary the SA from 1 to NSAs
%     array_x0 = arraystack_x (atemp, 1:Nx); %use atemp 10x 10
%     array_y0 = arraystack_y (atemp, 1:Ny); 
%  
%     SA_center= @(array_x)( min(array_x(:))+ ( max(array_x(:))-min(array_x(:)) )/2 );  
%     SA_x0 (ai) = SA_center (array_x0) ;
%     SA_y0 (ia) = SA_center (array_y0) ;     %y:oppositeOrder than x
%                  
%      ia = NSAs - ai;
%     atemp = Nx + atemp;         %10*ai + [1:1:Nx]; %next loop value
%  end
% end
% end
% %%
% [SA_x00 SA_y00] = meshgrid (SA_x0, SA_y0);
% figure(9)
% plot( SA_x00, SA_y00, 'vm'); grid on
% hold on
% plot(xall_SA, yall_SA, 'ob'); grid on


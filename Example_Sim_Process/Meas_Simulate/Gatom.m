function [gtot, g, g_shift] = Gatom(A,f,sigma,sfreq,t0,tend, shift_delay)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize Parameters
f_n = f*sfreq;                                                              %modulation
p = 1/sfreq;                                                                %sampling in time
sigma_n = sigma*p^2;                                                        %width
t = p:p:p*tend;                                                             %time vector
interp_st = 15;                                                             %interpolate last 15  pts of orig                            
interp_en = 8;                                                              %interpolate first 8 pts of shift

%% Generate Gabor Atom
g = A*exp(-1*(t-t0*p).^2/(2*sigma_n)).*cos(2*pi*f_n*(t-t0*p));

%% Generated Shifted Gabor Atom
if shift_delay > 0
    len_g = size(g,2);
    g_shift = zeros(1,len_g);
    shift_st = ceil(shift_delay/tend*len_g);
    shift_en = len_g;
    g_shift(:,shift_st:shift_en) =  g(:,1:shift_en-shift_st+1);             %no offset correction

    x = [shift_st-interp_st shift_st+interp_en];    v = [g_shift(:,shift_st-interp_st) g(:,interp_en+1)]; 
    xq = shift_st-interp_st:shift_st+interp_en;     g_shift_edge_smooth = interp1(x, v, xq, 'spline');
    g_shift(:,shift_st-interp_st:shift_st+interp_en) = smooth(g_shift_edge_smooth); 
else
    g_shift = zeros(size(g)); 
end

%% Total with Shift
gtot = g + g_shift;

%% Plot Gabor Atoms with Shift
%figure, set(gcf,'color','white')
%hold all; plot(gtot,'r'); plot(g,'k'); plot(g_shift,'g'); 
%legend('total','orig','shift')

end
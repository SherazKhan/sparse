function [sp, Us, Ss, Gstheta] = read_Gdecomp(sfwd, raw)
  if raw == 0
    sp = sfwd.numwhite_modes;                   %#white modes for 95% NMRA
    Us = [sfwd.whiteUs_p];                      %reduced singular vectors for white scort fwd
    Ss = [sfwd.whiteS_p];                       %reduced singular values for white scort fwd
    Gstheta = [sfwd.whiteGstheta_p];            %reduced rank approx to white scort fwd
    currstr = sfwd.currstr;                     %current strength that premultiplied G
    Gstheta = Gstheta/currstr;                  %return Gstheta without current strength
  else
    sp = sfwd.numraw_modes;                     %#raw modes
    Us = [sfwd.rawUs_p];                        %reduced singular vectors for raw scort fwd
    Ss = [sfwd.rawS_p];                         %reduced singular values for raw scort fwd
    Gstheta = [sfwd.rawGstheta_p];              %reduced rank approx to raw scort fwd   
    currstr = sfwd.currstr;                     %current strength that premultiplied G
    Gstheta = Gstheta/currstr;                  %return Gstheta without current strength
  end
end
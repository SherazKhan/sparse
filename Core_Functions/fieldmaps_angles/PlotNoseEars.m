function PlotNoseEars(hAxes, R, is2D)
    % Define coordinates
    NoseX = [0.8996; 1.15; 0.8996] * R/1.5;
    NoseY = [.14;       0;   -.14] * R/1.5;
    EarX  = [.0555 .0775 .0783 .0746  .0555  -.0055 -.0932 -.1313 -.1384 -.1199] * R * 1.5;
    EarY  = [.974, 1     1.016 1.0398 1.0638  1.06   1.074  1.044, 1      .958 ] * R;
    % 2D projection only
    if is2D
        scfactor = 2;
        EarX  = scfactor * EarX;
        EarY  = scfactor * EarY;
        NoseX = scfactor * NoseX;
        NoseY = scfactor * NoseY;
    end
        
    HLINEWIDTH = 2;
    HCOLOR = [.4 .4 .4];
    hold on
    % Plot nose
    plot(hAxes, -NoseY, 0.045+NoseX, ...
         'Color',     HCOLOR, ...
         'LineWidth', HLINEWIDTH, ...
         'tag',       'Nose');
    % Plot left ear
    plot(hAxes, -EarY, EarX, ...
         'Color',     HCOLOR, ...
         'LineWidth', HLINEWIDTH, ...
         'tag',       'leftEar');
    % Plot right ear
    plot(hAxes, EarY, EarX, ...
         'Color',     HCOLOR, ...
         'LineWidth', HLINEWIDTH, ...
         'tag',       'rightEar');
end
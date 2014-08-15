function PrintNearOptFigToPdf(hFigure,FileName)
    %Makes a pdf of the figure generated by nearoptplotmo.m with handle hFigure and saves to FileName
 
    %Set paper layout for the pdf print
    set(hFigure,'InvertHardcopy','off', 'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0.25 0.25 10.5 8], ...
            'PaperSize',[11 8.5]);
    print('-dpdf', '-r0', '-loose', FileName);
end


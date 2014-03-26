function saveFigure(figure, fileName, fileFormat)
    
    myStyle = hgexport('factorystyle');
    myStyle.Format = fileFormat;
   % myStyle.Width = 6;
    %myStyle.Height = 2.5;
    %myStyle.Resolution = 300;
    %myStyle.Units = 'inch';
    %myStyle.FixedFontSize = 12;
    hgexport(figure,fileName,myStyle,'Format',fileFormat)

end
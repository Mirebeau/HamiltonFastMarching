function [NewCoords] = RescaledCoords(OldCoords,Origin,Spacing)
    len = size(OldCoords); len=len(2);
    Ones = ones(1,len);
    Origins = Origin*Ones;
    Spacings = Spacing*Ones;
    
    NewCoords = 1+(OldCoords-Origins)./Spacings;
end
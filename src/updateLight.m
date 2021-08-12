function updateLight(src,event)
%UPDATELIGHT Summary of this function goes here
%   Detailed explanation goes here
l = findobj(gcf,'Type','light');
if isempty(l)
    camlight('headlight')
else
    for nr =1 : length(l)
        camlight(l(nr),'headlight')
    end
end

end


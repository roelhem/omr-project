function result = pickColors(colormap,n)
%PICKCOLORS Summary of this function goes here
%   Detailed explanation goes here
    picks = round(linspace(1,height(colormap), n));
    result = colormap(picks,:);
end


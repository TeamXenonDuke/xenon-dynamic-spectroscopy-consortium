function [] = boxplot_neat(all_data, g, colors, varargin)

% Create box plot with subjects marked with a dot '.'
boxplot(all_data,g,'Colors',colors,'Symbol','.','OutlierSize',10,'widths',0.8)

for nType = 1:max(g)
    subjs = find(g == nType);
    plot(g(subjs),all_data(subjs),'.','Color',colors(nType,:),'MarkerSize',14)
end 


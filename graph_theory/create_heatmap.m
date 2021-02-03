function [heatmap_fig] = create_heatmap(mat, xsize, ysize, xlab, ylab, format, title_str)

if ~exist('title_str','var')
     % title parameter does not exist, so default it to something
      title_str = '';
end
 
heatmap_fig = figure;
imagesc(mat);            % Create a colored plot of the matrix values
colormap(jet); 

textStrings = num2str(mat(:), format);       % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
[x, y] = meshgrid(1:xsize,1:ysize);  % Create x and y coordinates for the strings
hStrings = text(x(:), y(:), textStrings(:), ...  % Plot the strings
'HorizontalAlignment', 'center');
midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range

textColors = repmat(mat(:) > (midValue + (midValue/4)) | mat(:) < (midValue - (midValue/4)), 1, 3);  % Choose white or black for the
%   text color of the strings so
%   they can be easily seen over
%   the background color

set(hStrings, {'Color'}, num2cell(textColors, 2));  % Change the text colors
set(gca, 'XTick', 1:xsize)
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','FontName'),'FontName','Helvetica Neue')
title(title_str);
xlabel(xlab);
ylabel(ylab);
colorbar;

end
function figure_format(style)
%Formats the figure with some font and changes the size
%Dial 1 for full Page Width
%Dial 2 for half page width


if style==1
    den=1;
else
    den=2;
end
    

% Modify the text size and plot size
set(gca,'FontSize',14);
set(gcf,'Position',[100 100 1140/den 600])

end
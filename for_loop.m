figure                                  % Open a new figure                                   % End loop
for ii = 1:9                            % Start loop, have counter ii run from 1 to 9
    subplot(3,3,ii)                     % Draw into the subplot ii, arranged in 3 rows, 3 columns
    h = bar(1,1);                       % Fill the plot with a uniform color
    set(h,'FaceColor',[ii/9 ii/77 ii/1000]);      % Draw each in slightly different color
end                                     % End loop
% Many programs that are actually useful crucially depend on user input. This input
% comes mostly in one of two forms: from the mouse or from the keyboard.

figure                      % Opens a new figure
hold on;                    % Draw all lines in same figure
xlim([0 1])                 % Establish x-limits
ylim([0 1])                 % Y-limits
for ii = 1:200                % Start for-loop. Allow to draw 5 lines
    a = ginput(2);          % Get user input for 2 points
    plot(a(:,1),a(:,2));    % Draw the line.
end                         % The end.
%% Handle properties
%{ 
Handles are very useful. You can assign the property of an object (figure) to a
handle (ha). You can see their property through get(ha) and activate them
through set(ha)
%}
%1. Declare a handle
ha = figure;
%2. Check your variables
whos
%3. Get property of your handle (that you can control)
handle_property = get(ha)
%4. Activate your handle by calling their names
set(ha, 'Visible', 'off')
set(ha, 'Pointer', 'crosshair')
t = text(0.5, 0.5, input('Type something funny: ', 's'))
% Check out the property of your text (get) and assign (get)
text_property = get(t)
set(t, 'FontSize', 40, 'FontWeight', 'bold', 'Color', [0.4, 0.7, 0.8], 'HorizontalAlignment', 'center')
set(ha, 'Visible', 'on')        % I turned the visibiliy of the graph on again so that the complete product can be seen, rather than when its still in construction.
%5. The pause function allows user input to continue with single key
ii = 0;
disp('Try pressing a key: ')
% Pause MATLAB will stop all codes until user press a key
pause
% We can store the key input into a variable by get(x, 'CurrentCharacter')
ha_key = get(ha, 'CurrentCharacter')
% We can run it in a loop 9 times
while ii < 9
    pause
    if ii == 0
        disp('Press another key: ')
    end
    if ii >= 1
        disp('Continue to press until done...')
    end
    ha_key = get(ha, 'CurrentCharacter')
    ii = ii + 1;
end

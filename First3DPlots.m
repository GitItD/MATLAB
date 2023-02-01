a = -2:0.1:2;                  % Creating a vector with 21 elements (from -2.0 to 2.0 with step size of 0.2).
[x, y] = meshgrid(a, a*2);       % Creating x an y as a meshgrid of a
z = (sin(x)^2 + cos(y) - 1);         % Take the 2D exponential of x and y
figure                         % Open a new figure
MM = subplot(1,2,1)            % Create a left subplot
mesh(z)                        % Draw a wire mesh plot of data in z
TT = subplot(1,2,2)            % Create a right subplot 
surf_z = surf(z)                        % Draw a surface plot of data in z

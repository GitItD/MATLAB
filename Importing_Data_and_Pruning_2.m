%{ 
Data pruning
M1 = Data from movie, The Matrix 1; M2 = Data from the sequel, The Matrix 2 Reloaded. 
%}

%%

Movies = [];                                                % New Movies variable. Initializing it as empty.
temp = [M1, M2, M3];                                            % Create a joint temporary Matrix, consisting of two long vectors from imported data.
k = 1;                                                      % Initializing index k as 1
for ii = 1:length(temp)                                     % Could have said 1603 (the total numbers of rows in imported data)
    if isnan(temp(ii, 1)) == 0 & isnan(temp(ii, 2)) == 0 & isnan (temp(ii, 3)) == 0   % Only if both columns are NaN (not a number) is false (0), run the next 2 lines of code
        Movies(k,:) = temp(ii,:);                           % Fill the entries of 'temp' into new variable 'Movies'
        k = k + 1;                                          % Update k index only in this case
    end                                                     % End if clause
end                                                         % End for loop

Correlation_between_Matrix_1_2 = corrcoef(Movies(:,1),Movies(:,2))                           % Correlation p-values between M1, M2, and M3.
Correlation_between_Matrix_1_3 = corrcoef(Movies(:,1),Movies(:,3))
Correlation_between_Matrix_2_3 = corrcoef(Movies(:,2),Movies(:,3))
figure
hold on;
plot(Correlation_between_Matrix_1_2, 'color', 'k')
plot(Correlation_between_Matrix_1_3, 'color', 'g')
plot(Correlation_between_Matrix_2_3, 'color', 'r')
title('Correlations between the movies')
xlabel('Black = Corr Matrix 1 and 2, Green = Corr Matrix 1 and 3, Red = Corr Matrix 2 and 3')

figure                                                      % Create figure
hold on;
plot(Movies(:,1),Movies(:,2),'.','markersize',24)          % Plot ratings data of M1 vs M2 vs M3
plot(Movies(:,1),Movies(:,3),'.','markersize',24)
plot(Movies(:,2),Movies(:,3),'.','markersize',24)
title('Ratings for Matrix 1, 2, and 3')

figure                                                      % Open another new figure
subplot(3,1,1)                                              % Open new subplot
hold on;                                                    % Hold the plot
histogram(Movies(:,1),9)                                    % M1 data with 9 bins, which is enough since we only have 9 ratings
histfit(Movies(:,1),9)                                      % Fitting the histogram into a normalized Gaussian curve
xlim([-1 6])                                               % Create a plotting range
title('The Matrix Ratings')                                         % Add a title
subplot(3,1,2)                                              % Open a new subplot
hold on;                                                    % Hold the second subplot
histogram(Movies(:,2),9)                                    % Create histogram using M2 data using 9 bins
histfit(Movies(:,2),9)                                      % Fit M2 data into Gaussian
xlim([-1 6])                                                % Create plotting range
title('The Matrix Reloaded Ratings')                                % Add title
subplot(3,1,3)
histogram(Movies(:,3),9)
histfit(Movies(:,3),9)
xlim([-1 6])
title('The Matrix Revolution Ratings')

%%
%{
The following codes will improve the scatter plot by creating a 3D surface plot that provides more information:
%}

MT1 = (Movies(:,1)*2)+1;                                    % Assign a temp matrix in place of M1, multiply by 2 to get integral steps and add 1 matrix indices at the start w/1, not 0
MT2 = (Movies(:,2)*2)+1;                                    % Same thing but with M2.
MT3 = (Movies(:,3)*2)+1;
c = zeros(9,9);                                             % Create a matrix 'c' filled with 0 in 9X9
d = zeros(9,9);
e = zeros(9,9);
ii = 1;                                                     % Initialize index
for ii = 1:length(Movies)                                   % Start ii loop. This loop fills 'c' matrix with movie rating counts
    c(10-MT1(ii,1), MT2(ii,1)) = c(10-MT1(ii,1), MT2(ii,1)) + 1;     % Adding 1 in the cell count of 'c'
    d(10-MT1(ii,1), MT3(ii,1)) = d(10-MT1(ii,1), MT3(ii,1)) + 1;
    e(10-MT2(ii,1), MT3(ii,1)) = e(10-MT2(ii,1), MT3(ii,1)) + 1;
end                                                         % End loop

figure                                                      % New fig
surf(c)                                                     % Create surface of 'c'
shading interp                                              % Interpolate the shading
xlabel('Ratings for Matrix 1')                              % Label for the x-axis
ylabel('Ratings for Matrix 2')                              % Label for the y-axis
zlabel('Frequency')                                         % Label for the z-axis

figure
surf(d)
shading interp
xlabel('Ratings for Matrix 1')
ylabel('Ratings for Matrix 3')
zlabel('Frequency')

figure
surf(e)
shading interp
xlabel('Ratings for Matrix 2')
ylabel('Ratings for Matrix 3')
zlabel('Frequency')



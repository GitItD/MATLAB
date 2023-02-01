%% One-Way ANOVA example
hogg = ...
[24 14 11 7 19; 
    15 7 9 7 24; 
    21 12 7 4 19; 
    27 17 13 7 15; 
    33 14 12 12 10; 
    23 17 18 18 20];
[p, tbl, stats] = anova1(hogg)

stim_a = [39.2, 45.7, 45.9, 42.8, 60.2, 50.7, 39.9, 50.8, 43.0, 55.9];
stim_b = [43.2, 56.7, 32.8, 61.2, 54.6, 44.6, 53.2, 43.3, 35.1, 53.7];
stim_c = [66.5, 54.5, 62.6, 45.6, 46.8, 34.9, 53.3, 60.1, 69.7, 61.0];
stim_total = [stim_a; stim_b; stim_c]'
[p, tbl, stats] = anova1(stim_total);

%% Two-Way ANOVA example
car_mileage = ...
[33.3000   34.5000   37.4000; 
    33.4000   34.8000   36.8000; 
    32.9000   33.8000   37.6000; 
    32.6000   33.4000   36.6000; 
    32.5000   33.7000   37.0000; 
    33.0000   33.9000   36.7000];
nmbcars = 3; % Number of cars from each model
[p, tbl, stats] = anova2(car_mileage, nmbcars)
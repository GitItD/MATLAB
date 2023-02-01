P = [];                                        
a0 = 6000;
a1 = 6000;
m0 = 0;                                       
w0 = 300;                                     
P(1,:) = [m0, w0, a0, a1];
figure                                          
subplot(2,1,2), plot(0,a1,'.', 'markersize', 24);
title('Worker Numbers')
hold on;
xlim([0 25]) 
subplot(2,1,1), plot(0,a0,'.', 'markersize', 24);
title('Production Value')
hold on;
xlim([0 25])                                    
ii = 1;                                        
while ii < 25                                   
    P                                           
    a = input('How many workers this month? ')  
    b = 20 * a - a0
    c = a0*ii
    a0 = b;
    a1 = c;
    subplot(2,1,2), plot(ii, a1, '.', 'markersize', 24);
    subplot(2,1,1), plot(ii, a0, '.', 'markersize', 24);       
    P(ii + 1,:) = [ii, a, b, c];
    subplot(2,1,2), plot(P(:,1), P(:,4), 'color', 'k');
    title('Worker Numbers')
    subplot(2,1,1), plot(P(:,1), P(:,3), 'color', 'k');
    title('Production Value')
    ii = ii + 1;                                
end 

% Currently, we have both plots doing the same thing. How can I make it so
% that subplot 2,1,2 would display 1) dots of input "(a*ii)" and 2)
% connections of input "(a*ii)"?

% What you learned:
% By writing the problem in understandable format, you were able to create
% a program that solved the problem. 

% By thinking of an overall structure for how you would write the program,
% you were able to create an organized way of solving the problem, so
% whenever new problems arose, you were able to test to see where the new
% problem came from. In other words, you started the problem out with a
% simple assumption of what was correct, and then purposefully made
% mistakes, but only a single mistake each time, which made it easy to
% progress little by little until you found a way towards your goal by
% solving each new problem that arose on the way. 

% Whenever a new problem came, you came up with a hypothesis on why, and
% solved it by educated guesses. 

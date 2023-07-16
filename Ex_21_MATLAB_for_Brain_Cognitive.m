%% 21.7.4
%{
The theoretical scaling exponent of random white noise is 0.5. Can you
empirically confirm this, and does it depend on how much data you
have and how many scales you use? Run simulations for N data points
ranging from 1,000 to 1,000,000 and S scales ranging from 1% to 20%
of the length of the signal. Do you get a result close to 0.5? How about
if the data are not random noise, but a random series of +1 and ?1
(whatÅfs a good way to create a random series of +1 and ?1)?
%}

% My answer:

% Parameters
n = 1000;
wnoise = randn(1, n);
time = 1:length(n);
s_scales = ceil(logspace(log10(0.01*length(n)), log10(0.2*length(n)), n));
% Multiply FFT by scaling factor
fx = fft(wnoise);
fx = bsxfun(@times, fx, s_scales)';
y = real(ifft(fx));
% Calculate Hurst
y = y - mean(y);
Hurst = y./std(y);
% Plot Hurst over time
figure(1), clf
plot(time, Hurst);
xlabel('Time')
ylabel('Hurst Exponent')
title('White Noise Hurst')

% ChatGPT's answer (correction of my code):

% Parameters
n = 1000;
wnoise = randn(1, n);
time = 1:n; % Correction here
s_scales = ceil(logspace(log10(0.01*n), log10(0.2*n), n)); % Correction here
% Multiply FFT by scaling factor
fx = fft(wnoise);
fx = bsxfun(@times, fx, s_scales)';
y = real(ifft(fx));
% Calculate Hurst
y = y - mean(y);
Hurst = y./std(y);
% Calculate linear polynomial
p = polyfit(time', Hurst, 1);
yfit = polyval(p, time);
% Plot Hurst over time
figure(2), clf
hold on;
plot(time, Hurst);
plot(time, yfit, '--', 'color', 'r', 'lineWidth', 3);
xlabel('Time')
ylabel('Hurst Exponent')
title('ChatGPT Correction of Hurst Exponent of White Noise')

%% 21.7.7
%{
LetÅfs revisit the FWHM of the Gaussian in section 20.4 of chapter 20.
Can you compute this width more efficiently using some new tricks
you learned in this chapter?
%}
% We can compute this code below more efficiently by finding local maximum comparisons within 50% range from the mean instead.
% Create a Gaussian
x = -4:.1:4;
gauss = exp(-(x.^2));
gauss = gauss./max(gauss);
% Find maximum value and index
[max_val, max_idx] = max(gauss);
% Find half-maximum value
half_max = max_val/2;
% Find indices at which the Gaussian equals half-maximum value
left_idx = find(gauss(1:max_idx) < half_max, 1, 'last');
right_idx = find(gauss(max_idx:end) < half_max, 1) + max_idx - 1;
% Calculate FWHM
fwhm = x(right_idx) - x(left_idx);

%{
Original code:
% create a gaussian
x = -4:.1:4;
gaus = exp(-(x.^2));
gaus = gaus./max(gaus);
% find index of peak
[~,pidx] = max(gaus);
prepeak  = dsearchn(gaus(1:pidx)',.5);
postpeak = pidx-1+dsearchn(gaus(pidx:end)',.5);
% plot gaussian
figure(4), clf
plot(x, gaus, 'k-o')
hold on
% plot the empirical points closest to 50% amplitude
plot(x(prepeak),gaus(prepeak),'ro')
plot(x(postpeak),gaus(postpeak),'ro')
plot([ x(prepeak) x(postpeak)],[gaus(prepeak) gaus(postpeak)],'r--')
plot(get(gca,'xlim'),[.5 .5],'k:')
ylabel('Normalized amplitude')
%}
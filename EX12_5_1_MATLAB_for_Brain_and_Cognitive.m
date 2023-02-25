clear
clc

%% Ex 12.1, Convolve a signal in the time domain
signal = zeros(1,30);
signal(11:19) = 1;
kernel = [0, 0.31, 0.77, 1, 0.77, 0.31, 0];
kernel = kernel./sum(kernel);

N = length(signal);               % length time of the signal
M = length(kernel);             % length time of the kernel
halfKern = floor(M/2);       % Nyquist limit sampling rate, floor rounds M/2 toward negative infinity (ex: 1.9 = -2.000, 3.8 = 3)

dat4conv = [zeros(1, M-1), signal];
conv_res = zeros(1, N+M-1);

for ti = M+1:N-M+1
    tempdata = dat4conv(ti:ti+M-1);
    
    % compute dot-product (don't forget to flip the kernel backwards!)
    conv_res(ti) = sum( tempdata.*kernel(end:-1:1) );
end

w = conv(signal, kernel)

conv_res = conv_res(halfKern+1:end-halfKern);

figure(1), clf
hold on
plot(conv_res,'o-','linew',2,'markerface','g','markersize',9)
set(gca,'xlim',[0 N],'ylim',[min(signal)-.05 max(signal)+.05])
plot(w)
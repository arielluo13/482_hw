clear all; close all; clc;
%figure(1)
[y, Fs] = audioread('GNR.m4a');
y = y';
trgnr = length(y)/Fs; % record time in seconds
%plot((1:length(y))/Fs,y);
%xlabel('Time [sec]'); 
%ylabel('Amplitude');
%title('Sweet Child O'' Mine');
%p8 = audioplayer(y,Fs); playblocking(p8);
L = trgnr; n = length(y);
t2 = linspace(0,L,n+1); t = t2(1:n);
k = (1/L)*[0:(n/2 - 1) -n/2:-1]; ks = fftshift(k);
a = 50;
tau = 0:0.5:trgnr;
f = zeros(length(tau),1);
for j = 1:length(tau)
    g = exp(-a*(t - tau(j)).^2);
    yg = g.*y;
    ygt = fft(yg);
    ygtshift = abs(fftshift(ygt));
    x=ind2sub(size(ygtshift),find(ygtshift == max(ygtshift)));
    f(j)=ks(x(2));
end

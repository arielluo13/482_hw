clear all; close all; clc;
%figure(1)
[y, Fs] = audioread('Floyd.m4a');
y = y';
trgnr = length(y)/Fs; % record time in seconds
L = trgnr; n = length(y);
t2 = linspace(0,L,n+1); t = t2(1:n);
k = (1/L)*[0:(n/2) -n/2:-1]; ks = fftshift(k);
a = 50;
tau = 0:0.4:trgnr;
f = zeros(length(tau),1);

for j = 1:length(tau)
    g = exp(-a*(t - tau(j)).^2);
    yg = g.*y;
    ygt = fft(yg);
    ygtshift = abs(fftshift(ygt));
    x=ind2sub(size(ygtshift),find(ygtshift == max(ygtshift)));
    f(j)=ks(x(2));
end

for j = 1:length(f)
   if(f(j)>300)
       f(j)=0;
   end
end

filtered_f = zeros(length(tau),1);
for j = 1:length(tau)
    g = exp(-a*(t - tau(j)).^2);
    filter = exp(-0.2*((k - f(j)).^2));
    nyg = fft(y);
    nygt2 = filter.*nyg;
    inverse = abs(ifft(nygt2));
    inverse = g.*inverse;
    nygt = fft(inverse);
    nygt_spec(:,j)=fftshift(abs(nygt));
    nygtshift = abs(fftshift(nygt));
    ny=ind2sub(size(nygtshift),find(nygtshift == max(nygtshift)));
    filtered_f(j)=ks(ny);
end

  figure(1)
  pcolor(tau,ks,nygt_spec)
  shading interp
  set(gca,'ylim',[0,500],'Fontsize',16)
  colormap(hot)

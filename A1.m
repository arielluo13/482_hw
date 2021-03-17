%% Algorithm 1
% Clean workspace
clear all; close all; clc;

load subdata.mat % Imports the data as the 262144x49 (space by time) matrix called subdata 5

L = 10; % spatial domain
n = 64; % Fourier modes
x2 = linspace(-L,L,n+1); x = x2(1:n); y =x; z = x;
k = (2*pi/(2*L))*[0:(n/2 - 1) -n/2:-1]; ks = fftshift(k);

[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

Ut = zeros(n,n,n,49);
ave = zeros(n,n,n);
 for j=1:49
     Un(:,:,:)=reshape(subdata(:,j),n,n,n);
     ut = fftn(Un);
     Ut(:,:,:,j) = ut;
     ave = ave + ut;
 end
ave = abs(fftshift(ave))/49;
center = max(ave,[],'all')

%% Algorithm 2
[kx,ky,kz] = ind2sub(size(ave),find(ave == center));
k1 = Kx(kx,ky,kz); 
k2 = Ky(kx,ky,kz); 
k3 = Kz(kx,ky,kz);
[kx1,ky1,kz1]=meshgrid(k,k,k);
tau = 0.2;
filter = exp(-tau*((kx1 - k1).^2+(ky1 - k2).^2+(kz1 - k3).^2));
unft = zeros(64, 64, 64, 49);
unf = zeros(64, 64, 64, 49);
location = zeros(49,3);
for j = 1:49
    unft(:,:,:,j) = filter.*Ut(:,:,:,j);% apply filter
    inverse = abs(ifftn(unft(:,:,:,j)));
    unf(:,:,:,j) = inverse;
    peak = max(inverse,[],'all');
    [a1,b1,c1] = ind2sub(size(inverse),find(inverse == peak));
    a = X(a1,b1,c1); 
    b = Y(a1,b1,c1); 
    c = Z(a1,b1,c1);
    location(j,:) = [a,b,c];
    plot3(a,b,c,'.','Color','r')
    hold on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
%     pause(1)
end  
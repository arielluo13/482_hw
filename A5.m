vidObj=VideoReader('ski_drop_low.mp4')
mov=read(vidObj);

%% Code from lecture 27
% ut=fft(u);
% [t, utsol] = ode45(@(t,y) nls_rhs(t,y,k),t,ut);
% for j = 1:length(t)   
%     usol(j,:) = ifft(utsol(j,:)); % back to x-space
% end
% X = usol';
% X1 = X(:,1:end-1);
% X2 = X(:,2:end);
% [U, Sigma, V] = svd(X1,'econ');
% S = U'*X2*V*diag(1./diag(Sigma));
% [eV, D] = eig(S); % compute eigenvalues + eigenvectors
% mu = diag(D); % extract eigenvalues
% omega = log(mu)/dt;
% Phi = U*eV;

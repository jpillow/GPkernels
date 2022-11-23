% Examine periodic exponentiated-power-law kernel

addpath kernelfuns
addpath fouriertools

% 0. Set GP hyperparameters & grid points

rho = 2; % marginal variance
len = 5; % length scale
p = 1.5;  % power-law power

% Set x grid
nx = 40; % number of points
dx = 0.1;
xgrid = dx*(0:nx-1)'; % regular grid of x points
kperiod = nx*dx; % this is the period of the kernel


%% 1. Build periodic kernel explicitly in time domain

% % calculate # of periods needed to fall to 1e-8 of its height
% nperiods = min(100,ceil(((2*log(1e8))^(1/p))*len/kperiod));
% kexpowfun = @(x1,x2)(rho*exp(-0.5*abs((x1-x2)/len).^p));
% 
% K = zeros(nx,nx);
% % loop over periods
% for jj = -nperiods:nperiods
%     % Compute contribution to covariance 
%     K = K + kexpowfun(xgrid(:)+jj*kperiod,xgrid(:)'); % the covariance matrix
% end

% -------------------------
% EQUIVALENT FUNCTION CALL:
K = Kexppow_circular(xgrid,len,rho,p,kperiod);
% -------------------------

%% 2. Efficiently compute power spectrum and build in Fourier domain

% % efficiently compute Fourier spectrum
% jjperiod = -nperiods:nperiods;
% xgridAllLags = xgrid(:)+(jjperiod*kperiod);  % set of all lags (as a matrix)
% 
% % compute single row of K matrix
% krow = sum(rho*exp(-0.5*(abs(xgridAllLags)/len).^p),2); 
% 
% % compute its Fourier power spectrum
% cdiag = abs(fft(krow)); 
% 
% % Now make (inverse) FFT basis matrix
% Bfft = realfftbasis(nx)';
% 
% % now build K explicitly
% Kfd = Bfft*diag(cdiag)*Bfft';

% -------------------------
% EQUIVALENT FUNCTION CALL:
[cdiag,Bfft] = Kexppow_circular_fourierGridded(xgrid,len,rho,p,kperiod);
Kfd = Bfft*diag(cdiag2)*Bfft';
% -------------------------


%% 3.  Make plots

subplot(221); 
imagesc(K);  title(sprintf('time-domain K (p=%.2f)',p));

subplot(223);
imagesc(Kfd); title('Fourier-domain Kfd');

subplot(322); % plot some slices
iirow = [1,5,20]; % rows to plot
plot(xgrid, K(:,iirow),'-o',xgrid,Kfd(:,iirow),'--');
title('K slices');

% If desired, get the power spectrum of the original K using DFT matrix
cdiag_empir = diag(Bfft'*K*Bfft);

subplot(324); % plot difference between K and Kfd
plot(K-Kfd);  title('differences (K-Kfd)');

subplot(326);
plot(xgrid,cdiag_empir, '-o',xgrid,cdiag,'--x');
title('Eigenvalues (Fourier power spectrum)');
legend('from full K', 'from 1 row of K');
xlabel('cos / sin mode #');

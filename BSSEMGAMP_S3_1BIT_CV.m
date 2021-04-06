clc;  clear;   close all;
rng(6)      % 2,3,4,5,6 work
Monte = 500;
N = 181;             % grid size
maxit_outer = 1;
suppx = [ ]; DOAx = -90+((suppx-1)*180)/(N-1);
suppy = [ ]; DOAy = -90+((suppy-1)*180)/(N-1);
suppz = [ ]; DOAz = -90+((suppz-1)*180)/(N-1);
K = length(suppx);
x_dB = [20;20;20];  % amplitudes
theta_x = (-90:180/(N-1):90); theta_x_r = theta_x*pi/(N-1);
theta_y = (-90:180/(N-1):90); theta_y_r = theta_y*pi/(N-1);
theta_z = (-90:180/(N-1):90); theta_z_r = theta_z*pi/(N-1);
f_c = 1e9;  
w_lambda = 3e8 / f_c;  
d = 0.5 * w_lambda; % intersensor spacing
               

SNR = [];

SnrLen = length(SNR);
MSE_S3_1BIT_CV = zeros(1,SnrLen);
for SnrIdx = 1:SnrLen
    SNRdB = SNR(SnrIdx);
for monte = 1:Monte
L = 1;
c_sign = @(cpl_num)sign(real(cpl_num))+1j*sign(imag(cpl_num));

%% SLA
B1 = 101;B2 = 119;
q = [(0:B2-1)*B1, (1:2*B1-1)*B2];
xq= q.*d/w_lambda;
M = length(q);
M2 = max(q);

Ax = exp(-1j*2*pi*xq'*sin(theta_x_r)); % M*N
Ay = exp(-1j*2*pi*xq'*sin(theta_y_r)); % M*N
Az = exp(-1j*2*pi*xq'*sin(theta_z_r)); % M*N
x_amp = 10.^(x_dB/20);
x_amp = x_amp*ones(1,L);
rr = exp(1j*2*pi*rand(K,L));
Xx = zeros(N,L); Xx(suppx,:) = x_amp.* rr;%N*1
Xy = zeros(N,L); Xy(suppy,:) = x_amp.* rr;%N*1
Xz = zeros(N,L); Xz(suppz,:) = x_amp.* rr;%N*1
X = [Xx,Xy,Xz];
A = [Ax;Ay;Az];% 3M*N
if(L>1)
     wvar = ((norm(A*X,'fro'))^2/M/L)*10^(-SNRdB/10); 
else
     wvar = (norm(A*X))^2/M*10^(-SNRdB/10); 
end
% noise generation
w = sqrt(wvar/2)*randn(3*M,3*L)+1i*sqrt(wvar/2)*randn(3*M,3*L);
Y = c_sign(A*X+w);

tic
[Xhat, EMfin] = EMGMAMP(Y, A);
toc
t=toc;

[~,loc_peaks] = findpeaks( abs(Xhat(:,1)),'NPeaks',K,'SortStr','descend');
esti_theta_1 = sort(theta_x( loc_peaks ));
[~,loc_peaks] = findpeaks( abs(Xhat(:,2)),'NPeaks',K,'SortStr','descend');
esti_theta_2 = sort(theta_x( loc_peaks ));
[~,loc_peaks] = findpeaks( abs(Xhat(:,3)),'NPeaks',K,'SortStr','descend');
esti_theta_3 = sort(theta_x( loc_peaks ));
MSE_S3_1BIT_CV(SnrIdx) = MSE_S3_1BIT_CV(SnrIdx) + ...
   ((sum((esti_theta_1-DOAx).^2 + ...
      (esti_theta_2-DOAy).^2 + ...
      (esti_theta_3-DOAz).^2))/K)/3;


end
end
RMSE_S3_1BIT_CV = sqrt(MSE_S3_1BIT_CV/Monte);



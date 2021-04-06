clc;  clear;   close all;
rng(6)      % 2,3,4,5,6 work
Monte = 500;
N = 181;             % grid size
maxit_outer = 1;
suppx = [ ]; DOAx = -90+((suppx-1)*180)/(N-1);
suppy = [ ]; DOAy = -90+((suppy-1)*180)/(N-1);
suppz = [ ]; DOAz = -90+((suppz-1)*180)/(N-1);
fc1=2e3;      
fc2=3e3;
fc3=5e3;
K = length(suppx);
x_dB = [20;20;20];  % amplitudes
theta_x = (-90:180/(N-1):90); theta_x_r = theta_x*pi/180;
theta_y = (-90:180/(N-1):90); theta_y_r = theta_y*pi/180;
theta_z = (-90:180/(N-1):90); theta_z_r = theta_z*pi/180;
w_lambda = 3e8 / fc3;  % wavelength
d = 0.5 * w_lambda; % intersensor spacing               
SNR = [];

SnrLen = length(SNR);
MSE_S3_W_1BIT_C_CV = zeros(1,SnrLen);
for SnrIdx = 1:SnrLen
    SNRdB = SNR(SnrIdx);
for monte = 1:Monte
L = 1;
c_sign = @(cpl_num)sign(real(cpl_num))+1j*sign(imag(cpl_num));

%% SLA
B1 = 101;B2 = 119;
q = [(0:B2-1)*B1, (1:2*B1-1)*B2];
xq= q.*d;
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

A1 = exp(-1j*2*pi*fc1*xq'*sin(theta_x_r)); % M*N
A2 = exp(-1j*2*pi*fc2*xq'*sin(theta_x_r)); % M*N
A3 = exp(-1j*2*pi*fc3*xq'*sin(theta_x_r)); % M*N

wvar = 10^(-SNRdB/10); 
w = sqrt(wvar/2)*randn(M,L)+1i*sqrt(wvar/2)*randn(M,L);
A = [A1;A2;A3];
Yx = c_sign(A*Xx+w);
Yy = c_sign(A*Xx+w);
Yz = c_sign(A*Xx+w);

tic
[Xhatx, ~] = EMGMAMP(Yx, A);
[Xhaty, ~] = EMGMAMP(Yy, A);
[Xhatz, ~] = EMGMAMP(Yz, A);
toc
t=toc;

[~,loc_peaksx] = findpeaks( abs(Xhatx),'NPeaks',K,'SortStr','descend');esti_thetax = sort(theta_x( loc_peaksx ));
[~,loc_peaksy] = findpeaks( abs(Xhaty),'NPeaks',K,'SortStr','descend');esti_thetay = sort(theta_y( loc_peaksy ));
[~,loc_peaksz] = findpeaks( abs(Xhatz),'NPeaks',K,'SortStr','descend');esti_thetaz = sort(theta_z( loc_peaksz ));
MSE_S3_W_1BIT_C_CV(SnrIdx) = MSE_S3_W_1BIT_C_CV(SnrIdx) + ...
   ((sum((esti_thetax-DOAx).^2 + ...
      (esti_thetay-DOAy).^2 + ...
      (esti_thetaz-DOAz).^2))/K)/3;


end
end
RMSE_S3_W_1BIT_C_CV= sqrt(MSE_S3_W_1BIT_C_CV/Monte);

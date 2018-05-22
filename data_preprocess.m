% data preprocessing

%% fast fourier tranform
% load data
load('data_3intervals.mat');
load('theta3.mat');

t = p_interval3(:,1);
p = p_interval3(:,2);
Ts = (t(end)-t(1))/(length(t)-1);
Fs = 1/Ts;

% perform fast Fourier transform
nfftp = length(p);
nfft2p = 2.^nextpow2(nfftp);
fp = fft(p,nfftp);
fp = fp(1:nfft2p/2);
xfft = Fs .* (0: (nfft2p/2-1))/nfft2p;

% plot the frequency domain
plot(xfft, abs(fp/max(fp)));
xlim([0,1]);

% cut off frequency
cut_offp = 100/Fs/2;
orderp = 32;
hp = fir1(orderp, cut_offp);

% convolution
con = conv(delta_p, hp);

%% conventional filtering methods
% polynomial data spline
% reshape pressure vector to calculate grid nodes as mass averega point
matrix_p = reshape(p_interval(1:2000),80,25);
matrix_t = reshape(t(1:2000),80,25);
t_k = mean(matrix_t,1);
p_k = mean(matrix_p,1);
cs = spline(t_k,p_k);
p_t = ppval(cs,t);
delta_p = p-p_t;

% plot results from polynomial spline
yyaxis left
plot(t, delta_p); ylim([-0.3,0.3]); 
yyaxis right
plot(theta3(:,1)+0.05, theta3(:,2)); ylim([-1.1e-7, 1.1e-7]);
xlim([388,394]);

% Savitzky–Golay filter
% interval length chosen to be 159 in this case
delta_p = p-sgolayfilt(p,1,159);

% plot results from Savitzky–Golay filter
plot(t,delta_p); ylim([-0.3,0.3]); xlim([132,139]);
xlabel('Days since 7-7, 2003 0:00 UTC');
ylabel('\Delta p (kPa)');
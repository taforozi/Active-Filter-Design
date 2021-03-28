%% Band Elimination Filter - Inverse Chebyshev
% Aforozi Thomais
% AEM 9291

%% specifications
% aem = [a1,a2,a3,a4]
aem = [9,2,9,1];

f0 = 1.8*1000; %(kHz)
w0 = 2*pi*f0;
% w0 = sqrt(w1*w2) = sqrt(w3*w4) same

% Pass band (from 0 to ˘1 & from ˘2 to infinity)
f1 = 1200 + 25*(9-aem(4)); %Hz
w1 = 2*pi*f1;

f2 = f0^2/f1;
w2 = 2*pi*f2;

D = (1/1.8)*(f0^2 - f1^2)/f1;
f3 = (-D + sqrt(D^2 + 4*f0^2))/2;
w3 = 2*pi*f3;

f4 = f0^2/f3;
w4 = 2*pi*f4;

a_min = 30 - aem(3); %(dB)
a_max = 0.5 + aem(4)/18; %(dB)
Wp = 1;
Ws = (w2-w1)/(w4-w3);

% Bandwidth
bw = w2-w1;

%% Filter design
% filter order for normalized Wp = 1
n = acosh(sqrt((10^(a_min/10)-1)/(10^(a_max/10)-1)))/(acosh(Ws));
n = ceil(n); % n = 4

e = 1/(sqrt(10^(a_min/10)-1)); % (9-148)
alpha = (1/n)* asinh(1/e); %(9-149)

w_hp = 1/(cosh((1/n)*acosh(1/e))); % (9-139)

% Butterworth angles (deg) (n=4)
theta_1 = 22.5; % ch.9 p.13
theta_2 = -22.5;
theta_3 = 67.5;
theta_4 = -67.5;

% Chebyshev Poles
pole_1 = -sinh(alpha)*cosd(theta_1) + 1i*cosh(alpha)*sind(theta_1);
pole_2 = -sinh(alpha)*cosd(theta_2) + 1i*cosh(alpha)*sind(theta_2);
pole_3 = -sinh(alpha)*cosd(theta_3) + 1i*cosh(alpha)*sind(theta_3);
pole_4 = -sinh(alpha)*cosd(theta_4) + 1i*cosh(alpha)*sind(theta_4);

W12 = sqrt(real(pole_1)^2 + imag(pole_1)^2);
Q12 = 1/(2*cos(atan(imag(pole_1)/real(pole_1))));
% Q12 = W12/abs(2*real(pole_1)); same

W34 = sqrt(real(pole_3)^2 + imag(pole_3)^2);
Q34 = 1/(2*cos(atan(imag(pole_3)/real(pole_3))));
% Q34 = W34/abs(2*real(pole_3)); same

% Inverse Chebyshev Poles
inv_W12 = 1/W12;
inv_Q12 = Q12;

inv_W34 = 1/W34;
inv_Q34 = Q34;

% Pole Scaling 
inv_W12 = inv_W12 * Ws;
inv_W34 = inv_W34 * Ws;

% Zeros (9-130)
z1 = sec(pi/(2*n));
z2 = sec(3*pi/(2*n));

% Zero Scaling
Wz1 = z1 * Ws;
Wz2 = z2 * Ws;

% pole inversion ICH
s12 = 1/inv_W12;
s34 = 1/inv_W34;

% zero inversion 
inv_Wz1 = 1/Wz1;
inv_Wz2 = 1/Wz2;

% high pass poles
s12_hp = -s12/(2*Q12);
W12_hp = sqrt(s12^2 - s12_hp^2);

s34_hp = -s34/(2*Q34);
W34_hp = sqrt(s34^2 - s34_hp^2);

%% Pole tranformation - Geffe algorithm (ch.11 p.6)
%%%% pole pair s12_hp +- jW12_hp %%%%
S12_geffe = abs(s12_hp);
Wmega12_geffe = abs(W12_hp);
qc = w0/bw;

C12 = S12_geffe^2 + Wmega12_geffe^2;
D12 = 2*S12_geffe/qc;
E12 = 4 + C12/(qc^2);
G12 = sqrt(E12^2 - 4*D12^2);
Q12_geffe = (1/D12) * sqrt((1/2)*(E12+G12));
k12 = S12_geffe*Q12_geffe/qc;
W12_geffe = k12 + sqrt(k12^2 - 1);

w12_02 = W12_geffe*w0;
w12_01 = (1/W12_geffe)*w0;

% new poles
% 1.
p12_sigma1 = w12_01 / (2*Q12_geffe);
p12_wmega1 = w12_01 * sin(acos(1/(2*Q12_geffe)));

% 2.
p12_sigma2 = w12_02 / (2*Q12_geffe);
p12_wmega2 = w12_02 * sin(acos(1/(2*Q12_geffe)));

% +2 zeros at 0

%%%% pole pair -s34_hp +- jW34_hp %%%%
S34_geffe = abs(s34_hp);
Wmega34_geffe = abs(W34_hp);

C34 = S34_geffe^2 + Wmega34_geffe^2;
D34 = 2*S34_geffe/qc;
E34 = 4 + C34/(qc^2);
G34 = sqrt(E34^2 - 4*D34^2);
Q34_geffe = (1/D34) * sqrt((1/2)*(E34+G34));
k34 = S34_geffe*Q34_geffe/qc;
W34_geffe = k34 + sqrt(k34^2 - 1);

w34_02 = W34_geffe*w0;
w34_01 = (1/W34_geffe)*w0;

% new poles
% 1.
p34_sigma1 = w34_01 / (2*Q34_geffe);
p34_wmega1 = w34_01 * sin(acos(1/(2*Q34_geffe)));

% 2.
p34_sigma2 = w34_02 / (2*Q34_geffe);
p34_wmega2 = w34_02 * sin(acos(1/(2*Q34_geffe)));

% +2 zeros at 0

%% Zero transformation
% what arises from this process is two pairs of zeros 
% and two poles at zero

%%% inv_Wz1 %%%
K1 = 2 + inv_Wz1^2/(qc^2);
x1 = (K1 + sqrt(K1^2 - 4))/2;

wmega1_z1 = w0 * sqrt(x1);
wmega1_z2 = w0/sqrt(x1);

%%% inv_Wz2 %%%
K2 = 2 + inv_Wz2^2/(qc^2);
x2 = (K2 + sqrt(K2^2 - 4))/2;

wmega2_z1 = w0 * sqrt(x2);
wmega2_z2 = w0/sqrt(x2);

%% Zero-pole grouping and unit design
% a4 = 1 -> 7.21 & 7.23

%% Unit (…) (LPN) - (7.23)
% w12_01 < wmega1_z1

un1_w0 = w12_01;
un1_wz = wmega1_z1;
un1_Q = Q12_geffe;
un1_W0 = 1; un1_Wz = un1_wz/un1_w0;  

un1_C = 1/(2*un1_Q);
un1_R1 = 1;
un1_R4 = 1;
un1_R2 = 4*(un1_Q^2);
un1_R3 = (un1_Wz^2) / (2*un1_Q^2);
un1_R5 = (4*(un1_Q^2))/(un1_Wz^2 - 1);

un1_k_high = 1/(1 + un1_R3);
un1_k_low = un1_k_high*(un1_wz/un1_w0)^2;

un1_kf = un1_w0;
% klimakopoihsh gia C = 0.1ÏF (a3 = 9)
un1_km = un1_C/(0.1*10^(-6)*un1_kf);

un1_C_new = 0.1*10^(-6);
un1_R1_new = un1_km*un1_R1;
un1_R2_new = un1_km*un1_R2;
un1_R3_new = un1_km*un1_R3;
un1_R4_new = un1_km*un1_R4;
un1_R5_new = un1_km*un1_R5;

% sel 13/29
T1 = tf([un1_k_high 0 un1_k_high*un1_wz^2],[1 un1_w0/un1_Q un1_w0^2]);

%% Unit (……) (HPN) - (7.21)
% w12_02 > wmega1_z2

un2_w0 = w12_02;
un2_wz = wmega1_z2;
un2_Q = Q12_geffe;
un2_W0 = 1; un2_Wz = un2_wz/un2_w0;

un2_k1 = (un2_W0/un2_Wz)^2 - 1; % (7-135)
un2_k2 = ((2 + un2_k1)*un2_Q^2)/((2 + un2_k1)*un2_Q^2 + 1); %(7.136)

un2_k_high = un2_k2 * (un2_W0/un2_Wz)^2;

un2_R1 = 1;
un2_R2 = un2_Q^2 *(un2_k1 + 2)^2;
un2_R3 = 1;
un2_R4 = un2_Q^2 *(un2_k1 + 2);
un2_C = 1/(un2_Q*(2 + un2_k1)); 
un2_C1 = un2_k1*un2_C;

un2_kf = un2_w0;
% klimakopoihsh gia C = 0.1ÏF (a3 = 9)
un2_km = un2_C/(un2_kf*0.1*10^(-6));

un2_C_new = 0.1*10^(-6);
un2_R1_new = un2_km * un2_R1;
un2_R2_new = un2_km * un2_R2;
un2_R3_new = un2_km * un2_R3;
un2_R4_new = un2_km * un2_R4;
un2_C1_new = un2_k1*un2_C_new;


% (7-134) ﬁ
T2 = tf([un2_k_high 0 un2_k_high*un2_wz^2],[1 un2_w0/un2_Q un2_w0^2]);


%% Unit (………) (LPN) - (7.23)
% w34_01 < wmega2_z1

un3_w0 = w34_01;
un3_wz = wmega2_z1;
un3_Q = Q34_geffe;
un3_W0 = 1; un3_Wz = un3_wz/un3_w0;  

un3_C = 1/(2*un3_Q);
un3_R1 = 1;
un3_R4 = 1;
un3_R2 = 4*un3_Q^2;
un3_R3 = un3_Wz^2 / (2*un3_Q^2);
un3_R5 = (4*un3_Q^2)/(un3_Wz^2 - 1);

un3_k_high = 1/(1 + un3_R3);
un3_k_low = un3_k_high*(un3_wz/un3_w0)^2;

un3_kf = un3_w0;
% scaling for C = 0.1ÏF (a3 = 9)
un3_km = un3_C/(0.1*10^(-6)*un3_kf);

un3_C_new = 0.1*10^(-6);
un3_R1_new = un3_km*un3_R1;
un3_R2_new = un3_km*un3_R2;
un3_R3_new = un3_km*un3_R3;
un3_R4_new = un3_km*un3_R4;
un3_R5_new = un3_km*un3_R5;

T3 = tf([un3_k_high 0 un3_k_high*un3_wz^2],[1 un3_w0/un3_Q un3_w0^2]);

%% Unit (…V) (HPN) - (7.21)
% w34_02 > wmega2_z2

un4_w0 = w34_02;
un4_wz = wmega2_z2;
un4_Q = Q34_geffe;
un4_W0 = 1; un4_Wz = un4_wz/un4_w0;

un4_k1 = (un4_W0/un4_Wz)^2 - 1; % (7-135)
un4_k2 = ((2 + un4_k1)*un4_Q^2)/((2 + un4_k1)*un4_Q^2 + 1); % (7-136)

un4_k_high = un4_k2 * (un4_W0/un4_Wz)^2; % (1-137)

un4_R1 = 1;
un4_R2 = un4_Q^2 *(un4_k1 + 2)^2;
un4_R3 = 1;
un4_R4 = un4_Q^2 *(un4_k1 + 2); % (7-139)
un4_C = 1/(un4_Q*(2 + un4_k1)); % (7-140)
un4_C1 = un4_k1*un4_C;

un4_kf = un4_w0;
% scaling for C = 0.1ÏF (a3 = 9)
un4_km = un4_C/(un4_kf*0.1*10^(-6));

un4_C_new = 0.1*10^(-6);
un4_R1_new = un4_km * un4_R1;
un4_R2_new = un4_km * un4_R2;
un4_R3_new = un4_km * un4_R3;
un4_R4_new = un4_km * un4_R4;
un4_C1_new = un4_k1*un4_C_new;

T4 = tf([un4_k_high 0 un4_k_high*un4_wz^2],[1 un4_w0/un4_Q un4_w0^2]);

%% Gain setting
% desired low frequency gain 10 dB (a4 = 1)

total_gain = un1_k_high*un2_k_high*un3_k_high*un4_k_high;
% 20log(k) = 10 dB -> k = 10^(10/20)

k_des = 10^(10/20); % > total_gain -> gain amplification is needed
a = k_des / total_gain; % > 1 
% 1 + r2/r1 = a and r1 = 10kOhm
% r2 = (a - 1)*r1 = 4394 Ohm = 4.394 kOhm;

T = T1*T2*T3*T4;
Ta = a*T;

%% Plots
ltiview({'bodemag'}, Ta, T1, T2, T3, T4);

plot_transfer_function(T1,[f0 f1 f2 f3 f4]);
plot_transfer_function(T2,[f0 f1 f2 f3 f4]);
plot_transfer_function(T3,[f0 f1 f2 f3 f4]);
plot_transfer_function(T4,[f0 f1 f2 f3 f4]);
plot_transfer_function(Ta,[f0 f1 f2 f3 f4]);

% ¡¸Û‚ÂÛÁ
plot_transfer_function(inv(Ta), [f0 f1 f2 f3 f4]);

plot_transfer_function(inv(T), [10 f0 f1 f2 f3 f4]);

%% Periodic input signal
% a4 = 1 -> ·)
T_period = 20*1/2000; 
f_deigm = 10^6; % sampling frequency
dt_deigm = 1/f_deigm;
t_end = T_period - dt_deigm;
time = 0:dt_deigm:t_end;

u = 0.8*cos((w0 - (w0 - w3)/2)*time) + 1.0*cos((w0 + (w0 + w3)/2)*time) ...
    + cos(0.5*w1*time) + 0.8*cos(2.4*w2*time) + 0.4*cos(3.5*w2*time);

% Frequencies for ac sources in multisim
fu_1 = (w0 - (w0 - w3)/2)/(2*pi); % 0.8Vpk
fu_2 = (w0 + (w0 + w3)/2)/(2*pi); % 1.0Vpk
fu_3 = (0.5*w1)/(2*pi); % 1.0Vpk
fu_4 = 2.4*w2/(2*pi); % 0.8Vpk
fu_5 = 3.5*w2/(2*pi); % 0.4Vpk

figure;
plot(time,u,'LineWidth',0.8);
grid on;
axis([0 inf -4 5]);
title('Input periodic signal','Interpreter','Latex');
xlabel('Time (sec)','Interpreter','Latex');
ylabel('Amplitude','Interpreter','Latex');


y = lsim(Ta,u,time);
figure;
plot(time,y,'LineWidth',0.8,'Color','red');
grid on;
axis([0 inf -10 10]);
title('Output signal','Interpreter','Latex');
xlabel('Time (sec)','Interpreter','Latex');
ylabel('Amplitude','Interpreter','Latex');

figure;
plot(time,u,time,y,'LineWidth',0.8);
legend('u','y')
grid on;
axis([0 0.007 -10 10]);
title('Input-Output signals','Interpreter','Latex');
xlabel('Time (sec)','Interpreter','Latex');
ylabel('Amplitude','Interpreter','Latex');
%% Fourier

U_fft = fft(u);
Lu = length(u);

P2_U = abs(U_fft/Lu);
P1_U = P2_U(1:Lu/2+1);
P1_U(2:end-1) = 2*P1_U(2:end-1);

fu = f_deigm*(0:(Lu/2))/Lu;
figure;
plot(fu, P1_U);
xticks([0 750 2000 3500 5000 10000 15000]);
axis([0 1.5*10^4 0 1]);
title('Single-Sided Amplitude Spectrum of U(t)','Interpreter','Latex');
xlabel('f (Hz)','Interpreter','Latex');
ylabel('P1 $$u$$ (f)','Interpreter','Latex');


Y_fft = fft(y);
Ly = length(y);

P2_Y = abs(Y_fft/Ly);
P1_Y = P2_Y(1:Ly/2+1);
P1_Y(2:end-1) = 2*P1_Y(2:end-1);

fy = f_deigm*(0:(Ly/2))/Ly;
figure;
plot(fy, P1_Y);
xticks([0 500 750 1500 2000 2500 3000 3500 4000 4500 5000 10000 15000]);
axis([0 5000 0 3.5]);
title('Single-Sided Amplitude Spectrum of Y(t)','Interpreter','Latex');
xlabel('f (Hz)','Interpreter','Latex');
ylabel('P1 $$y$$ (f)','Interpreter','Latex');


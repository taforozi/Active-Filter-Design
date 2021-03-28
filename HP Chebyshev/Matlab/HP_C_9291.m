%% High Pass Filter - Chebyshev
% Aforozi Thomais
% AEM 9291

%% specifications
% aem = [a1,a2,a3,a4];

aem = [9,2,9,1];
m = 1; % a3 = 9;

% Pass Band Frequency (kHz)
fp = (3 + m)*1000; 
wp = 2*pi*fp;

% Stop Band Frequency (kHz)
fs = fp/1.8;
ws = 2*pi*fs;

% Damping
a_min = 25 + aem(3) * (4/9); % (dB)
a_max = 0.5 + aem(4) * (0.25/9); % (dB)

Wp = 1;
Ws = wp/ws;

%% Filter design
% Order
n = acosh(sqrt((10^(a_min/10)-1)/(10^(a_max/10)-1)))/(acosh(Ws));
n = ceil(n);

% Half power frequency
e = sqrt(10^(a_max/10) - 1); % (9-76)
alpha = (1/n)* asinh(1/e); %(9-149)
%  Whp = cosh((1/n)*acosh(1/sqrt(10^(a_max/10)-1)));
Whp = cosh((1/n)*acosh(1/e));

% poles & zeros
% Butterworth angles (deg) (n=5)
theta_1 = 0;
theta_2 = 36;
theta_3 = -36;
theta_4 = 72;
theta_5 = -72;

% Chebyshev poles
pole_1 = -sinh(alpha)*cosd(theta_1) + 1i*cosh(alpha)*sind(theta_1);
pole_2 = -sinh(alpha)*cosd(theta_2) + 1i*cosh(alpha)*sind(theta_2);
pole_3 = -sinh(alpha)*cosd(theta_3) + 1i*cosh(alpha)*sind(theta_3);
pole_4 = -sinh(alpha)*cosd(theta_4) + 1i*cosh(alpha)*sind(theta_4);
pole_5 = -sinh(alpha)*cosd(theta_5) + 1i*cosh(alpha)*sind(theta_5);

% Ω, Q
W0_1 = sqrt(real(pole_1)^2 + imag(pole_1)^2); 
W0_23 = sqrt(real(pole_2)^2 + imag(pole_2)^2);
W0_45 = sqrt(real(pole_4)^2 + imag(pole_4)^2);

% Q1 = W0_1/abs(2*real(pole_1)); same
% Q23 = W0_23/abs(2*real(pole_2));
% Q45 = W0_45/abs(2*real(pole_4));

Q1 =  1/(2*cos(atan(imag(pole_1)/real(pole_1))));
Q23 =  1/(2*cos(atan(imag(pole_2)/real(pole_2))));
Q45 =  1/(2*cos(atan(imag(pole_4)/real(pole_4))));


% Pole inversion -> Highpass poles

whp = wp/Whp;
f_hp = whp/(2*pi);

s1 = wp/abs(pole_1); % (rad/s) ω0 unit (Ι) Q1
w23 = wp/W0_23; % (rad/s) ω0 unit (ΙΙ) Q23
w45 = wp/W0_45; % (rad/s) ω0 unit (ΙΙΙ) Q45

%% Unit design
%% Unit (Ι) - Simple CR circuit
% T1(s) = s/(s+p1)

un1_w0 = s1;
un1_kf = un1_w0;

% scaling for C = 0.1μF = 0.1*10^(-6) = 10^(-7) (a4 = 1)
un1_C = 1;
un1_R1 = 1;
un1_km = 1/(un1_kf*10^(-7));

un1_R1_new = un1_km * un1_R1;
un1_C_new = 10^(-7);

T1 = tf([1 0],[1 1/(un1_R1_new*un1_C_new)]);

%% Unit (ΙΙ) - Highpass Sallen Key Strategy (2) (α3=9)

un2_w0 = w23;
un2_Q = Q23;

un2_C1 = 1; % (6-76)
un2_C2 = 1;
un2_k = 1;
un2_R1 = 1/(2*un2_Q); % (6-77)
un2_R2 = 2*un2_Q;

un2_kf = un2_w0;

% scaling for C = 0.1*10^(-6)
un2_km = 1/(un2_kf*10^(-7));

un2_C1_new = 10^(-7);
un2_C2_new = 10^(-7);
un2_R1_new = un2_km * un2_R1;
un2_R2_new = un2_km * un2_R2;

% (6-60)
T2 = tf([un2_k 0 0],[1 un2_w0/un2_Q un2_w0^2]);

%% Μονάδα (ΙΙΙ) - Highpass Sallen Key Strategy(2) 

un3_w0 = w45;
un3_Q = Q45;

un3_C1 = 1; % (6-76)
un3_C2 = 1;
un3_k = 1;
un3_R1 = 1/(2*un3_Q); % (6-77)
un3_R2 = 2*un3_Q;

un3_kf = un3_w0;

% scaling for C = 0.1*10^(-6)
un3_km = 1/(un3_kf*10^(-7));

un3_C1_new = 10^(-7);
un3_C2_new = 10^(-7);
un3_R1_new = un3_km * un3_R1;
un3_R2_new = un3_km * un3_R2;

% (6-60)
T3 = tf([un3_k 0 0],[1 un3_w0/un3_Q un3_w0^2]);

%% Gain setting - desired high frequency gain 10dB
% total high frequency gain
total_gain = 1*un2_k*un3_k;

% desired 10dB -> 20*log(x) = 10 -> x = 10^(0.5)
a = 10^(0.5)/total_gain; % -> gain amplification is needed
% 1 + r2/r1 = a kai epilegw r1 = 10kOhm
% r2 = (a - 1)*r1 = 21.6228 kOhm;

T = T1*T2*T3;
Ta = a*T;

%% plots

plot_transfer_function(T1,[fs f_hp fp]);
plot_transfer_function(T2,[fs f_hp fp]);
plot_transfer_function(T3,[fs f_hp fp]);
% 10^6 (υψηλές συχνότητες για να φανεί το κέρδος)
plot_transfer_function(T,[fs f_hp fp 10^6]);
plot_transfer_function(Ta,[fs f_hp fp 10^6]);

ltiview({'bodemag'}, Ta, T1, T2, T3);

plot_transfer_function(inv(T),[fs f_hp fp 10^6]);
plot_transfer_function(inv(Ta),[fs f_hp fp 10^6]);


%% Periodic input signal
% α4 = 1 -> β)
T_period = 20*1/2000; 
f_deigm = 10^6; % sampling frequency
dt_deigm = 1/f_deigm;
t_end = T_period - dt_deigm;
time = 0:dt_deigm:t_end;

u = cos(0.4*ws*time) + 0.5*cos(0.9*ws*time) ...
    + cos(1.4*wp*time) + 0.7*cos(2.4*wp*time) + 0.5*cos(4.5*wp*time);

% Frequencies for ac sources in multisim
fu_1 = (0.4*ws)/(2*pi); % 1Vpk
fu_2 = (0.9*ws)/(2*pi); % 0.5Vpk
fu_3 = 1.4*wp/(2*pi); % 1Vpk
fu_4 = 2.4*wp/(2*pi); % 0.7Vpk
fu_5 = 4.5*wp/(2*pi); % 0.5Vpk


figure;
plot(time,u,'LineWidth',0.8);
grid on;
axis([0 0.005 -4 4]);
title('Input periodic signal','Interpreter','Latex');
xlabel('Time (sec)','Interpreter','Latex');
ylabel('Amplitude','Interpreter','Latex');


y = lsim(Ta,u,time);
figure;
plot(time,y,'LineWidth',0.8,'Color','red');
grid on;
axis([0 0.005 -8 12]);
title('Output signal','Interpreter','Latex');
xlabel('Time (sec)','Interpreter','Latex');
ylabel('Amplitude','Interpreter','Latex');

figure;
plot(time,u,time,y,'LineWidth',0.8);
legend('u','y')
grid on;
axis([0 0.005 -10 15]);
title('Input-Output signals','Interpreter','Latex');
xlabel('Time (sec)','Interpreter','Latex');
ylabel('Amplitude','Interpreter','Latex');

%% Input/Output Spectrum

U_fft = fft(u);
Lu = length(u);

P2_U = abs(U_fft/Lu);
P1_U = P2_U(1:Lu/2+1);
P1_U(2:end-1) = 2*P1_U(2:end-1);

fu = f_deigm*(0:(Lu/2))/Lu;
figure;
plot(fu, P1_U);
axis([0 4*10^4 0 inf]);
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
axis([0 25000 0 inf]);
title('Single-Sided Amplitude Spectrum of Y(t)','Interpreter','Latex');
xlabel('f (Hz)','Interpreter','Latex');
ylabel('P1 $$y$$ (f)','Interpreter','Latex');



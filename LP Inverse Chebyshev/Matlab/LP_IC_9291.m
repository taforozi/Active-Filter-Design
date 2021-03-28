%% LowPass Filter - Inverse Chebyshev

%% specifications
% aem = [a1,a2,a3,a4]
aem = [9,2,9,1];
m = 1; 

% Pass Band Frequency (kHz)
fp = 1.1*(3+m)*1000;
wp = 2*pi*fp;

% Stop Band Frequency (kHz)
fs = 1.7*fp;
ws = 2*pi*fs;

% Damping (dB)
a_min = 25 + (max(1,aem(3))-5)*(3/4);
a_max = 0.5 + (max(1,aem(4))-5)/16;

% Normalization -> ws = 1
Ws = 1; 
Wp = wp/ws;

%% Filter design
% 1. Filter order (9-137) for normalized Ws = 1 
% n = acosh(((10^(a_min/10)-1)/(10^(a_max/10)-1))^(1/2))/acosh(1/Wp);  
n = acosh(sqrt((10^(a_min/10)-1)/(10^(a_max/10)-1)))/acosh(1/Wp); % (9-137)
n = ceil(n); % n = 5 

% 2. Half power frequency
e = 1/(sqrt(10^(a_min/10)-1)); % (9-148), (9-123)
% e = sqrt(10^(a_max/10) - 1); % (9-76)
w_hp = 1/(cosh((1/n)*acosh(1/e))); % (9-139)
% f_hp = w_hp/(2*pi);
f_hp = w_hp * fp; % w_hp *wp/(2*pi)

alpha = (1/n)* asinh(1/e); %(9-149)

% Butterworth angles (deg) n=5
theta_1 = 0;
theta_2 = 36;
theta_3 = -36;
theta_4 = 72;
theta_5 = -72;

%% Poles & zeros
% chebyshev (9-102),(9-103)
% since we have angles in deg we use sind & cosd
pole_1 = -sinh(alpha)*cosd(theta_1) + 1i*cosh(alpha)*sind(theta_1);
pole_2 = -sinh(alpha)*cosd(theta_2) + 1i*cosh(alpha)*sind(theta_2);
pole_3 = -sinh(alpha)*cosd(theta_3) + 1i*cosh(alpha)*sind(theta_3);
pole_4 = -sinh(alpha)*cosd(theta_4) + 1i*cosh(alpha)*sind(theta_4);
pole_5 = -sinh(alpha)*cosd(theta_5) + 1i*cosh(alpha)*sind(theta_5);

% inv chebyshev
% angle(pole) = (rad) -> cos,sin
inv_p1 = 1/(abs(pole_1))*(cos(angle(pole_1))+ 1i*sin(angle(pole_1)));
inv_p2 = 1/(abs(pole_2))*(cos(angle(pole_2))+ 1i*sin(angle(pole_2)));
inv_p3 = 1/(abs(pole_3))*(cos(angle(pole_3))+ 1i*sin(angle(pole_3)));
inv_p4 = 1/(abs(pole_4))*(cos(angle(pole_4))+ 1i*sin(angle(pole_4)));
inv_p5 = 1/(abs(pole_5))*(cos(angle(pole_5))+ 1i*sin(angle(pole_5)));

w01 = sqrt(real(inv_p1)^2 + imag(inv_p1)^2);
w02 = sqrt(real(inv_p2)^2 + imag(inv_p2)^2); % gia tous polous 2 3
% w03 = sqrt(real(inv_p3)^2 + imag(inv_p3)^2)
w04 = sqrt(real(inv_p4)^2 + imag(inv_p4)^2); % gia tous polous 4 5
% sqrt(real(inv_p5)^2 + imag(inv_p5)^2);

Q1 = 1/(2*cos(atan(imag(inv_p1)/real(inv_p1))));
Q2 = 1/(2*cos(atan(imag(inv_p2)/real(inv_p2))));
Q4 = 1/(2*cos(atan(imag(inv_p4)/real(inv_p4))));

% same
% Q1 = w01/abs(2*real(inv_p1));
% Q2 = w02/abs(2*real(inv_p2));
% Q4 = w04/abs(2*real(inv_p4));

% zeros (9-143)
z1 = sec(pi/(2*n));
z2 = sec(3*pi/(2*n));
% z3 = sec(5*pi/(2*n)); -> apeiro

%% Poles - zeros grouping and LPN filter design 
% a3 = 9 -> 7.23
%% Unit (É) (LPN) - poles 4,5 - zero z1 
un1_w0 = w04;
un1_wz = z1; % wz>w0
un1_Q = Q4;
un1_W0 = 1; un1_Wz = un1_wz/un1_w0; % >1

un1_C = 1/(2*un1_Q); % (7-150)
un1_R2 = 4*(un1_Q)^2; % (7-150)
un1_R5 = (4*(un1_Q)^2)/((un1_Wz)^2 -1); % (7-152)
un1_R3 = (un1_Wz)^2/(2*(un1_Q)^2); % (7-155)
un1_R1 = 1; un1_R4 = 1;

k1_high = 1/(1+un1_R3); % (7-154)
k1_low = k1_high*(un1_wz/un1_w0)^2;

un1_kf = ws*un1_w0;

% scaling for C = 0.1ìF
un1_km = un1_C/(ws*un1_w0*10^(-7));

un1_C_new = 10^(-7);
un1_R1_new = un1_km * un1_R1;
un1_R2_new = un1_km * un1_R2;
un1_R3_new = un1_km * un1_R3;
un1_R4_new = un1_km * un1_R4;
un1_R5_new = un1_km * un1_R5;

% frequency scaling
un1_wz = un1_wz*ws;
un1_w0 = un1_w0*ws;
% (7-121)
T1 = tf([k1_high 0 k1_high*un1_wz^2],[1 (un1_w0/un1_Q) un1_w0^2]);

%% Unit (ÉI) (LPN) - poles 2,3 - zero z2 (in a similar way) 
un2_w0 = w02;
un2_wz = z2; % wz>w0
un2_Q = Q2;
un2_W0 = 1; un2_Wz = un2_wz/un2_w0; % >1

un2_C = 1/(2*un2_Q); % (7-150)
un2_R2 = 4*(un2_Q)^2; % (7-150)
un2_R5 = (4*(un2_Q)^2)/((un2_Wz)^2 -1); % (7-152)
un2_R3 = (un2_Wz)^2/(2*(un2_Q)^2); % (7-155)
un2_R1 = 1; un2_R4 = 1;

k2_high = 1/(1+un2_R3); % (7-154)
k2_low = k2_high*(un2_wz/un2_w0)^2;

un2_kf = ws*un2_w0;

% scaling for C = 0.1ìF
un2_km = un2_C/(ws*un2_w0*10^(-7));

un2_C_new = 10^(-7);
un2_R1_new = un2_km * un2_R1;
un2_R2_new = un2_km * un2_R2;
un2_R3_new = un2_km * un2_R3;
un2_R4_new = un2_km * un2_R4;
un2_R5_new = un2_km * un2_R5;

% frequency scaling
un2_wz = un2_wz*ws;
un2_w0 = un2_w0*ws;
% (7-121)
T2 = tf([k2_high 0 k2_high*un2_wz^2],[1 (un2_w0/un2_Q) un2_w0^2]);

%% Unit (ÉÉÉ) (first order) - Real pole

un3_w0 = w01;
un3_C = 1;
un3_R = 1/un3_w0;

un3_kf = ws;
un3_km = un3_C/(10^(-7)*un3_kf);

un3_C_new = 10^(-7);
un3_R_new = un3_km*un3_R;

T3 = tf(un3_w0*ws,[1 un3_w0*ws]); %(2-9)

%% Gain setting
% a4 = 1 -> desired gain = 0dB or 1
k = k1_low*k2_low*1;
k_des = 1; % 0 dB
% 20log(aK) = 0
gain = k_des/k;

% gain = -r2/r1 = 0.7673
% if r1 = 10kOhm -> r2 = 7.673 kOhm
r1 = 10000;
r2 = gain*r1;

T = gain*T1*T2*T3;

T_before = T1 * T2 * T3;
invTbefore = inv(T_before);

%% Plots
% Response
ltiview({'bodemag'}, T, T1, T2, T3);

plot_transfer_function(T1, [fp f_hp fs]);
plot_transfer_function(T2, [fp f_hp fs]);
plot_transfer_function(T3, [fp f_hp fs]);
plot_transfer_function(T_before, [10 fp f_hp fs]);

plot_transfer_function(T, [10 fp f_hp fs]);

% Damping
plot_transfer_function(inv(T), [10 fp f_hp fs]);

plot_transfer_function(invTbefore, [10 fp f_hp fs]);

%% Periodic input signal 
% a4 = 1 -> pulse wave ô = 0.4Ô
% basic frequency 2kHz

f_sig = 2*1000;

T_period = 10* 1/f_sig;  
f_deigm = 10^6; % sampling frequency
dt_deigm = 1/f_deigm;
t_end = T_period - dt_deigm;
time = 0:dt_deigm:t_end;

% input signal
% duty cycle 40%
u = 1*(square(2*pi*2000*(time),40)+1)/2;  

figure;
plot(time,u,'LineWidth',0.8);
grid on;
axis([0 inf -0.5 1.5]);
title('Input periodic signal','Interpreter','Latex');
xlabel('Time (sec)','Interpreter','Latex');
ylabel('Amplitude','Interpreter','Latex');


y = lsim(T,u,time);
figure;
plot(time,y,'LineWidth',0.8,'Color','red');
grid on;
axis([0 inf -1 1.5]);
title('Output signal','Interpreter','Latex');
xlabel('Time (sec)','Interpreter','Latex');
ylabel('Amplitude','Interpreter','Latex');


figure;
plot(time,u,time,y,'LineWidth',0.8);
legend('u','y')
grid on;
axis([0 0.005 -0.5 2]);
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
axis([0 80000 0 inf]);
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
axis([0 10^4 0 inf]);
title('Single-Sided Amplitude Spectrum of Y(t)','Interpreter','Latex');
xlabel('f (Hz)','Interpreter','Latex');
ylabel('P1 $$y$$ (f)','Interpreter','Latex');



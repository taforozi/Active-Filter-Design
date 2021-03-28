%% BandPass Filter - Butterworth


%% specifications
% aem = [a1,a2,a3,a4]
aem = [9,2,9,1];

f0 = 1.3 * 1000; %(kHz)
w0 = 2*pi*f0; % (rad/sec) 
% w0 = sqrt(w1*w2) same

% Pass Band
f1 = 750 + 25 * aem(4); %(Hz)
w1 = 2*pi*f1;

f2 = (f0^2)/f1;
w2 = 2*pi*f2;

D = 2.7 * (f0^2 - f1^2)/f1;

bw = w2-w1;

% Stop Band
f3 = (-D + sqrt(D^2 + 4 * f0^2))/2;
w3 = 2*pi*f3;

f4 = (f0^2)/f3;
w4 = 2*pi*f4;

a_min = 25 + aem(3) * (5/9); % (dB)
a_max = 0.5 + aem(4)/36; % (dB)

Wp = 1;
Ws = (w4-w3)/(w2-w1);


%% Filter design
% Filter order (9-52) for Wp = 1
n = (log((10^(a_min/10)-1)/(10^(a_max/10)-1)))/(2*log(Ws/Wp));
n = ceil(n);

% (9-48) - Frequency at 3dB 
W0 = Wp / ((10^(a_max/10)-1)^(1/(2*n))); %(rad/sec)

% Butterworth angles (deg) n=5
theta_1 = 0;
theta_2 = 36;
theta_3 = -36;
theta_4 = 72;
theta_5 = -72;

%% Butterworth poles
s1 = W0 * (-cosd(theta_1) + 1i*sind(theta_1));
s2 = W0 * (-cosd(theta_2) + 1i*sind(theta_2));
s3 = W0 * (-cosd(theta_3) + 1i*sind(theta_3));
s4 = W0 * (-cosd(theta_4) + 1i*sind(theta_4));
s5 = W0 * (-cosd(theta_5) + 1i*sind(theta_5));

%% Pole transformation (Geffe algorithm)
%%%%%% s1: real pole %%%%%%
qc = w0/bw;
S1 = abs(s1);

Q1 = qc/S1;
psi_1 = acos(1/(2*Q1)); % angle in rad
psi_1 = rad2deg(psi_1); % angle in deg

% 2 poles arise from this process which lie on the cycle
% with radius = w0 
p1 = w0 * (-cosd(psi_1) + 1i*sind(psi_1));
p2 = w0 * (-cosd(psi_1) - 1i*sind(psi_1));

%%%%%% s2,3: first complex pair of poles %%%%%
S2 = abs(real(s2)); 
Wmega2 = abs(imag(s2));

C23 = S2^2 + Wmega2^2;
D23 = (2*S2)/qc;
E23 = 4 + C23/(qc^2);
G23 = sqrt(E23^2 - 4*D23^2);
Q23 = (1/D23)*sqrt((1/2)*(E23 + G23));
k23 = (S2*Q23)/qc;
W23 = k23 + sqrt(k23^2 - 1);

w23_02 = W23*w0; % (11-35) 
w23_01 = (1/W23)*w0; 

% angles psi
psi_23 = acosd(1/(2*Q23));

% first pair of poles
s1_23 = w23_01/(2*Q23);
w1_23 = w23_01*sind(psi_23);

% second pair of poles
s2_23 = w23_02/(2*Q23);
w2_23 = w23_02*sind(psi_23);

%%%%%% s4,5: second complex pair of poles %%%%%%
S4 = abs(real(s4));
Wmega4 = abs(imag(s4));

C45 = S4^2 + Wmega4^2;
D45 = (2*S4)/qc;
E45 = 4 + C45/(qc^2);
G45 = sqrt(E45^2 - 4*D45^2);
Q45 = (1/D45)*sqrt((1/2)*(E45 + G45));
k45 = (S4*Q45)/qc;
W45 = k45 + sqrt(k45^2 - 1);

w45_02 = W45*w0; 
w45_01 = (1/W45)*w0;

psi_45 = acosd(1/(2*Q45));

% first pair of poles
s1_45 = w45_01/(2*Q45);
w1_45 = w45_01*sind(psi_45);

% second pair of poles
s2_45 = w45_02/(2*Q45);
w2_45 = w45_02*sind(psi_45);

%% Units: a3 = 9 -> Strategy (2) (Delyiannis-Fried)
% kef 7 sel 22-23
% Unit I: poles of s1 -> w0,Q1
% Unit II: poles of s2,3 -> w1,Q23
% Unit III: poles of s2,3 -> w2,Q23
% Unit IV: poles of s4,5 -> w1,Q45
% Unit V: poles of s4,5 -> w2,Q45

%% Unit (I) 
un1_Q = Q1;
un1_w0 = w0;
un1_R1 = 1/(2*un1_Q);
un1_R2 = 2*un1_Q;
un1_C1 = 1;
un1_C2 = 1;

% scaling for C = 1.0ìF = 10^(-6) (a2 = 2)
un1_kf = w0;
un1_km = un1_C1/(w0*10^(-6));

un1_R1_new = un1_km * un1_R1;
un1_R2_new = un1_km * un1_R2;
un1_C1_new = un1_C1/(w0*un1_km);
un1_C2_new = un1_C1_new;

T1 = tf([(-2*un1_Q*un1_w0) 0], [1 (un1_w0/un1_Q) un1_w0^2]); 

%% Unit (II)
un2_Q = Q23;
un2_w0 = w23_01;
un2_kf = un2_w0;

un2_R1 = 1/(2*un2_Q);
un2_R2 = 2*un2_Q;
un2_C1 = 1;
un2_C2 = 1;

un2_km = 1/(un2_w0*10^(-6));

un2_C1_new = 10^(-6);
un2_C2_new = 10^(-6);
un2_R1_new = un2_R1 * un2_km; 
un2_R2_new = un2_R2 * un2_km; 

% (11-57)
T2 = tf([(-2*un2_Q*un2_w0) 0], [1 (un2_w0/un2_Q) un2_w0^2]); 

%% Unit (III)
un3_Q = Q23; 
un3_w0 = w23_02;
un3_kf = un3_w0;

un3_C1 = 1;
un3_C2 = 1;
un3_R1 = 1/(2*un3_Q);
un3_R2 = 2*un3_Q;

un3_km = 1/(un3_w0*10^(-6));

un3_C1_new = 1/(un3_km*un3_kf); %C/(kf*km) = 10^(-6)
un3_C2_new = 1/(un3_km*un3_kf);
un3_R1_new = un3_km*un3_R1;
un3_R2_new = un3_km*un3_R2;

% (11-57)
T3 = tf([(-2*un3_Q*un3_w0) 0], [1 (un3_w0/un3_Q) un3_w0^2]); 

%% Unit (IV)
un4_Q = Q45;
un4_w0 = w45_01;
un4_kf = un4_w0;

un4_C1 = 1;
un4_C2 = 1;
un4_R1 = 1/(2*un4_Q);
un4_R2 = 2*un4_Q;

un4_km = 1/(un4_w0*10^(-6));

un4_C1_new = 1/(un4_km*un4_kf);
un4_C2_new = 1/(un4_km*un4_kf);
un4_R1_new = un4_km*un4_R1;
un4_R2_new = un4_km*un4_R2;

% (11-57)
T4 = tf([(-2*un4_Q*un4_w0) 0], [1 (un4_w0/un4_Q) un4_w0^2]); 

%% Unit (V)
un5_Q = Q45;
un5_w0 = w45_02;
un5_kf = un5_w0;

un5_C1 = 1;
un5_C2 = 1;
un5_R1 = 1/(2*un5_Q);
un5_R2 = 2*un5_Q;

un5_km = 1/(un5_w0*10^(-6));

un5_C1_new = 1/(un5_km*un5_kf);
un5_C2_new = 1/(un4_km*un4_kf);
un5_R1_new = un5_km*un5_R1;
un5_R2_new = un5_km*un5_R2;

% (11-57)
T5 = tf([(-2*un5_Q*un5_w0) 0], [1 (un5_w0/un5_Q) un5_w0^2]); 

%% Gain setting - Desired gain 10dB (a4 = 1)
% total response gain 
% |T1(jù0)|*|T2(jù0)|*|T3(jù0)|*|T4(jù0)|*|T5(jù0)|
% (11-58) or norm(evalfr(Tf, i*ù0)) for each tf
T1_norm = sqrt((2*un1_Q*un1_w0*w0)^2/((un1_w0^2-w0^2)^2 +((un1_w0/un1_Q)^2)*w0^2));
T2_norm = sqrt((2*un2_Q*un2_w0*w0)^2/((un2_w0^2-w0^2)^2 +((un2_w0/un2_Q)^2)*w0^2));
T3_norm = sqrt((2*un3_Q*un3_w0*w0)^2/((un3_w0^2-w0^2)^2 +((un3_w0/un3_Q)^2)*w0^2));
T4_norm = sqrt((2*un4_Q*un4_w0*w0)^2/((un4_w0^2-w0^2)^2 +((un4_w0/un4_Q)^2)*w0^2));
T5_norm = sqrt((2*un5_Q*un5_w0*w0)^2/((un5_w0^2-w0^2)^2 +((un5_w0/un5_Q)^2)*w0^2));

total_gain = T1_norm*T2_norm*T3_norm*T4_norm*T5_norm;
% Caution! The result total_gain isn't in dB
a = (10^(10/20))/total_gain;

T = a*T1*T2*T3*T4*T5;

% resistors for gain setting
Z2 = un1_R1_new/a; Z3 = un1_R1_new/(1-a);

%% Plots
plot_transfer_function(T1, [f0 f1 f2 f3 f4]);
plot_transfer_function(T2, [f0 f1 f2 f3 f4]);
plot_transfer_function(T3, [f0 f1 f2 f3 f4]);
plot_transfer_function(T4, [f0 f1 f2 f3 f4]);
plot_transfer_function(T5, [f0 f1 f2 f3 f4]);
plot_transfer_function(T, [f0 f1 f2 f3 f4]);

ltiview({'bodemag'}, T1, T2, T3, T4, T5,T);

invT = inv(T);
plot_transfer_function(invT, [f0 f1 f2 f3 f4]);

T_before = T1*T2*T3*T4*T5;
plot_transfer_function(T_before, [f0 f1 f2 f3 f4]);
inv_Tb = inv(T_before);
plot_transfer_function(inv_Tb, [f0 f1 f2 f3 f4]);

%% Periodic input signal - Fourier Spectrum
% a4 = 1 -> ã)

T_period = 20*1/2000; 
f_deigm = 10^6; % sampling frequency
dt_deigm = 1/f_deigm;
t_end = T_period - dt_deigm;
time = 0:dt_deigm:t_end;

u = cos((w0 - (w0-w1)/2)*time) + 0.8*cos((w0 + (w0+w1)/3)*time) ...
    + 0.8*cos(0.4*w3*time) + 0.6*cos(2.5*w4*time) + 0.5*cos(3*w4*time);

% Frequencies for ac sources in multisim
fu_1 = (w0 - (w0-w1)/2)/(2*pi); % 1Vpk
fu_2 = (w0 + (w0+w1)/3)/(2*pi); % 0.8Vpk
fu_3 = 0.4*w3/(2*pi); % 0.8Vpk
fu_4 = 2.5*w4/(2*pi); % 0.6Vpk
fu_5 = 3*w4/(2*pi); % 0.5Vpk


figure;
plot(time,u,'LineWidth',0.8);
grid on;
axis([0 inf -4 4]);
title('Input periodic signal','Interpreter','Latex');
xlabel('Time (sec)','Interpreter','Latex');
ylabel('Amplitude','Interpreter','Latex');


y = lsim(T,u,time);
figure;
plot(time,y,'LineWidth',0.8,'Color','red');
grid on;
axis([0 inf -6 6]);
title('Output signal','Interpreter','Latex');
xlabel('Time (sec)','Interpreter','Latex');
ylabel('Amplitude','Interpreter','Latex');


figure;
plot(time,u,time,y,'LineWidth',0.8);
legend('u','y')
grid on;
axis([0 0.008 -6 6]);
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
axis([0 1.5*10^4 0 0.9]);
xticks([0 2000 5000 10000 15000]);
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
axis([0 5000 0 2.5]); 
title('Single-Sided Amplitude Spectrum of Y(t)','Interpreter','Latex');
xlabel('f (Hz)','Interpreter','Latex');
ylabel('P1 $$y$$ (f)','Interpreter','Latex');


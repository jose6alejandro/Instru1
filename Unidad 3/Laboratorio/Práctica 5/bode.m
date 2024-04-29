% Practica 5 - Bode
% Castro Jose, Hurtado Carlos

clc
clear all

global RF RL R1 R2 R3 C L

RF = 50;
R1 = 150 + RF; 
RL = 89;  
R2 = 1200 + RL;
R3 = 300;       
C = 1e-6;
L = 141.0e-3;

% Funciones de transferencia
syms s
ft_VC = (L*s+R2)/((C*L*R1+C*L*R3)*s^2+(C*R1*R2+C*R1*R3+C*R2*R3+L)*s+R1+R2);
ft_IL = (C*R3*s+1)/((C*L*R1+C*L*R3)*s^2+(C*R1*R2+C*R1*R3+C*R2*R3+L)*s+R1+R2);

% Funciones complejas
% VC
num = [L, R2];
dem = [(C*L*R1+C*L*R3), (C*R1*R2+C*R1*R3+C*R2*R3+L), R1+R2];
Gs_VC = tf(num, dem)

% IL
num = [C*R3, 1];
dem = [(C*L*R1+C*L*R3), (C*R1*R2+C*R1*R3+C*R2*R3+L), R1+R2];
Gs_IL = tf(num, dem);

% BODE PARA VC

% Matlab bode
[mag,phase,w] = bode(Gs_VC); 
magdB = 20*log10(squeeze(mag)); % Convertir la magnitud a decibelios
phaseDeg= squeeze(phase);

% Maple magnitud VC
magnitud_vc_maple = 20*log(sqrt(1661521+0.19881e-1*w.^2)./sqrt((-0.70500e-4*w.^2+1489).^2+.7148702500*w.^2))./log(10);

% Maple fase VC
fase_vc = ((.141*1i)*w+1289)./(-0.70500e-4*w.^2+(.845500*1i)*w+1489);
fase_vc_maple = atan(imag(fase_vc)./real(fase_vc)) * 57.2957795;

% Graficas
figure(1)
subplot(2,1,1);
semilogx(w, magdB,'g', w, magnitud_vc_maple, 'r');
title('Respuesta en Frecuencia de VC(s) - Magnitud'); 
xlabel('Frecuencia (rad/seg)'); 
ylabel('Magnitud (dB)'); 
legend('Matlab','Maple')
grid on; 

subplot(2,1,2);
semilogx(w, phaseDeg,'g', w, fase_vc_maple, 'r');
title('Respuesta en Frecuencia de VC(s) - Fase');
xlabel('Frecuencia (rad/seg)'); 
ylabel('Fase (deg)'); 
legend('Matlab','Maple') 
grid on; 

% BODE PARA IL

% Matlab bode
[mag,phase,w] = bode(Gs_IL); 
magdB = 20*log10(squeeze(mag)); % Convertir la magnitud a decibelios
phaseDeg= squeeze(phase);

% Maple magnitud IL
magnitud_il_maple = 20*log(sqrt(1+9.0000*10^(-8)*w.^2)./sqrt((-0.70500e-4*w.^2+1489).^2+.7148702500*w.^2))./log(10);

% Maple fase IL
fase_IL = ((0.300e-3*1i)*w+1)./(-0.70500e-4*w.^2+(.845500*1i)*w+1489);
fase_IL_maple = atan(imag(fase_IL)./real(fase_IL)) * 57.2957795;

% Graficas
figure(2)
subplot(2,1,1);
semilogx(w, magdB,'g', w, magnitud_il_maple, 'r');
title('Respuesta en Frecuencia de IL(s) - Magnitud'); 
xlabel('Frecuencia (rad/seg)'); 
ylabel('Magnitud (dB)'); 
legend('Matlab','Maple')
grid on; 

subplot(2,1,2);
semilogx(w, phaseDeg,'g', w, fase_IL_maple, 'r');
title('Respuesta en Frecuencia de IL(s) - Fase');
xlabel('Frecuencia (rad/seg)'); 
ylabel('Fase (deg)'); 
legend('Matlab','Maple') 
grid on; 

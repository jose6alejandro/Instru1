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
Gs_VC = tf(num, dem);

% IL
num = [C*R3, 1];
dem = [(C*L*R1+C*L*R3), (C*R1*R2+C*R1*R3+C*R2*R3+L), R1+R2];
Gs_IL = tf(num, dem);

% BODE PARA VC

% Medidas
magnitud = [
-1.110346557
-1.681455766
-2.18040806
-3.87640052
-5.352124804
-5.81460078
-8.13427866
-9.473214452
-11.83520069
-12.04119983
-14.70364354
-15.08974664
-15.91760035
-16.24958558
-16.83275016
-17.58852138
];

frecuencia = [
80
161
200
304
406
502
666
809
1100
1310
1504
1724
1907
2218
2601
2806
];

fase = [
-14.112
-26.4
-30.24
-43.63636364
-50.60240964
-52.46231156
-59.61038961
-67.31707317
-71.44432194
-73.098
-75.8016
-74.4768
-75.5172
-77.45256
-72.09972
-67.68072
];

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
semilogx(w/(2*pi), magdB,'g', w/(2*pi), magnitud_vc_maple, '--r', frecuencia, magnitud, 'b');
title('Respuesta en Frecuencia de VC(s) - Magnitud'); 
xlabel('Frecuencia (Hz)'); 
ylabel('Magnitud (dB)'); 
legend('Matlab','Maple', 'Experimental')
grid on; 

subplot(2,1,2);
semilogx(w/(2*pi), phaseDeg,'g', w/(2*pi), fase_vc_maple, '--r', frecuencia, fase, 'b');
title('Respuesta en Frecuencia de VC(s) - Fase');
xlabel('Frecuencia (Hz)'); 
ylabel('Fase (deg)'); 
legend('Matlab','Maple', 'Experimental') 
grid on; 

% BODE PARA IL

% Medidas
frecuencia = [
108
200
302
417
523
646
737
848
922
1120
1321
1511
1721
2013
2235
2400
2959
3294
];

magnitud = [
-63.8234965
-64.5887865
-65.30309755
-66.05435779
-66.64203334
-67.36627679
-67.62176255
-68.00271104
-68.43666191
-69.02576149
-69.65774122
-70.16385852
-70.51843386
-71.4738895
-71.88758147
-72.10205878
-73.51368027
-74.31602473
];

fase = [
-9.72
-16.056
-21.744
-28.5228
-30.1248
-31.3956
-37.1448
-38.77056
-38.1708
-40.7232
-45.1782
-48.9564
-49.5648
-56.88738
-61.1496
-63.936
-63.9144
-61.66368
];

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
semilogx(w/(2*pi), magdB,'g', w/(2*pi), magnitud_il_maple, '--r', frecuencia, magnitud, 'b');
title('Respuesta en Frecuencia de IL(s) - Magnitud'); 
xlabel('Frecuencia (Hz)'); 
ylabel('Magnitud (dB)'); 
legend('Matlab','Maple', 'Experimental')
grid on; 

subplot(2,1,2);
semilogx(w/(2*pi), phaseDeg,'g', w/(2*pi), fase_IL_maple, 'r', frecuencia, fase, 'b');
title('Respuesta en Frecuencia de IL(s) - Fase');
xlabel('Frecuencia (Hz)'); 
ylabel('Fase (deg)'); 
legend('Matlab','Maple', 'Experimental') 
grid on; 

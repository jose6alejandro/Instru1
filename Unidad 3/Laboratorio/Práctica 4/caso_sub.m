% Practica 4 - Segundo Orden
% Castro Jose, Hurtado Carlos

clc
clear all

global V R1 R2 R3 C L FREQ

V = 4.9;
RF = 50;
R1 = 153 + RF; 
RL = 89;  
R2 = 50 + RL;
R3 = 303;       
C = 1e-6;
L = 141.0e-3;

FREQ = 102.1; 
tiempo_simulacion = (1/FREQ)/2;
tiempo_total = 0;
tspan = [0 tiempo_simulacion];
ti = [0];

% Condiciones Iniciales
vC_inicial    = -R2*V/(R1+R2);
iL_inicial    = -V/(R1+R2); 
iC_inicial    = 2*V/(R1+R3);
vL_inicial    = 2*V*R3/(R1+R3);

vC_inicial_diff = 2*V/(C*(R1+R3));
iC_inicial_diff = -2*V*(C*R1*R3+L)/(C*(R1+R3)^2*L);
iL_inicial_diff = 2*V*R3/((R1+R3)*L);
vL_inicial_diff  = -2*V*(C*R1*R2*R3+C*R1*R3^2+C*R2*R3^2-L*R1)/(C*(R1+R3)^2*L);

condiciones_iniciales = [vC_inicial, iL_inicial, iC_inicial, vL_inicial];
valores_ode_total = condiciones_iniciales;

VC_maple_total = vC_inicial;
IC_maple_total = iC_inicial;
IL_maple_total = iL_inicial;
VL_maple_total = vL_inicial;

for i=1:6
    
    options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-4 1e-4]);
    [tiempo_ode, valores_ode] = ode45(@ecuaciones, tspan, condiciones_iniciales, options);
    
    valores_ode_total = [valores_ode_total; valores_ode];
    tiempo_total = [tiempo_total; tiempo_ode];
     
    
    %Maple
    tiempo_maple = (tiempo_ode - ti); 
    
    VC_maple = condiciones_iniciales(1).*cosh((1/2)*tiempo_maple*sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2)/(C*(R1+R3)*L)).*exp(-(1/2)*(C*R1*R2+C*R1*R3+C*R2*R3+L)*tiempo_maple/(C*(R1+R3)*L))+(-cosh((1/2)*tiempo_maple*sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2)/(C*(R1+R3)*L)).*exp(-(1/2)*(C*R1*R2+C*R1*R3+C*R2*R3+L)*tiempo_maple/(C*(R1+R3)*L))+1)*R2*V/(R1+R2)+sinh((1/2)*tiempo_maple*sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2)/(C*(R1+R3)*L)).*exp(-(1/2)*(C*R1*R2+C*R1*R3+C*R2*R3+L)*tiempo_maple/(C*(R1+R3)*L))*(2*C*L*vC_inicial_diff*R1*(R1+R2)+2*C*L*vC_inicial_diff*R3*(R1+R2)+C*condiciones_iniciales(1)*R1*R2*(R1+R2)+C*condiciones_iniciales(1)*R1*R3*(R1+R2)+C*condiciones_iniciales(1)*R2*R3*(R1+R2)-C*R1*R2^2*V-C*R2*R3*V*(R1+R2)+L*condiciones_iniciales(1)*(R1+R2)-L*R2*V)/(sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2)*(R1+R2));
    IC_maple = exp(-(1/2)*(C*R1*R2+C*R1*R3+C*R2*R3+L)*tiempo_maple/(C*(R1+R3)*L)).*(condiciones_iniciales(3).*cosh((1/2)*tiempo_maple*sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2)/(C*(R1+R3)*L))+sinh((1/2)*tiempo_maple*sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2)/(C*(R1+R3)*L))*(2*iC_inicial_diff*C*L*R1+2*iC_inicial_diff*C*L*R3+C*condiciones_iniciales(3)*R1*R2+C*condiciones_iniciales(3)*R1*R3+C*condiciones_iniciales(3)*R2*R3+L*condiciones_iniciales(3))/sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2));
    VL_maple = exp(-(1/2)*(C*R1*R2+C*R1*R3+C*R2*R3+L)*tiempo_maple/(C*(R1+R3)*L)).*(condiciones_iniciales(4).*cosh((1/2)*tiempo_maple*sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2)/(C*(R1+R3)*L))+sinh((1/2)*tiempo_maple*sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2)/(C*(R1+R3)*L))*(2*vL_inicial_diff*C*L*R1+2*vL_inicial_diff*C*L*R3+C*condiciones_iniciales(4)*R1*R2+C*condiciones_iniciales(4)*R1*R3+C*condiciones_iniciales(4)*R2*R3+L*condiciones_iniciales(4))/sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2));
    IL_maple = condiciones_iniciales(2).*cosh((1/2)*tiempo_maple*sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2)/(C*(R1+R3)*L)).*exp(-(1/2)*(C*R1*R2+C*R1*R3+C*R2*R3+L)*tiempo_maple/(C*(R1+R3)*L))+(-cosh((1/2)*tiempo_maple*sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2)/(C*(R1+R3)*L)).*exp(-(1/2)*(C*R1*R2+C*R1*R3+C*R2*R3+L)*tiempo_maple/(C*(R1+R3)*L))+1)*V/(R1+R2)+sinh((1/2)*tiempo_maple*sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2)/(C*(R1+R3)*L)).*exp(-(1/2)*(C*R1*R2+C*R1*R3+C*R2*R3+L)*tiempo_maple/(C*(R1+R3)*L))*(2*C*L*iL_inicial_diff*R1*(R1+R2)+2*C*L*iL_inicial_diff*R3*(R1+R2)+C*condiciones_iniciales(2)*R1*R2*(R1+R2)+C*condiciones_iniciales(2)*R1*R3*(R1+R2)+C*condiciones_iniciales(2)*R2*R3*(R1+R2)-C*R1*R2*V-C*R3*V*(R1+R2)+L*condiciones_iniciales(2)*(R1+R2)-L*V)/(sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2)*(R1+R2));
    
    VC_maple_total = [VC_maple_total; VC_maple];
    IC_maple_total = [IC_maple_total; IC_maple];
    VL_maple_total = [VL_maple_total; VL_maple];
    IL_maple_total = [IL_maple_total; IL_maple];
    
    ti = tiempo_total(end);
    V = -V;
    
    vC_inicial_diff = 2*V/(C*(R1+R3));
    iC_inicial_diff = -2*V*(C*R1*R3+L)/(C*(R1+R3)^2*L);
    iL_inicial_diff = 2*V*R3/((R1+R3)*L);
    vL_inicial_diff  = -2*V*(C*R1*R2*R3+C*R1*R3^2+C*R2*R3^2-L*R1)/(C*(R1+R3)^2*L);
    
    tspan = [tiempo_ode(end) tiempo_ode(end) + tiempo_simulacion];
        
    condiciones_iniciales = [valores_ode(end, 1) valores_ode(end, 2) valores_ode(end, 1)/(-R1-R3)+R1*valores_ode(end, 2)/(-R1-R3)-V/(-R1-R3) valores_ode(end, 1)+(valores_ode(end, 1)/(-R1-R3)+R1*valores_ode(end, 2)/(-R1-R3)-V/(-R1-R3))*R3-valores_ode(end, 2)*R2];
end

vc = importdata('Medidas/sub_capacitor/VOLTAJE.CSV',',',18);
voltaje_Data =  vc.data(:, 4:5);

ic = importdata('Medidas/sub_capacitor/CORRIENTE.CSV',',',18);
corriente_Data =  ic.data(:, 4:5);

il = importdata('Medidas/sub_inductor/CORRIENTE.CSV',',',18);
corriente2_Data =  il.data(:, 4:5);

vl = importdata('Medidas/sub_inductor/VOLTAJE.CSV',',',18);
voltaje2_Data =  vl.data(:, 4:5);

% Capacitor
figure(1);
subplot(2, 2, 1);
plot(tiempo_total,valores_ode_total(:,1),'g', tiempo_total, VC_maple_total,'r', voltaje_Data(:, 1) + 0.01465, voltaje_Data(:, 2),'b');
title('Voltaje del capacitor (sobreamortiguado)')
xlabel('Tiempo (seg)')
ylabel('Voltaje (V)')
legend('ODE','Maple','Experimental')
grid on;

subplot(2, 2, 2);
plot(tiempo_total, valores_ode_total(:,3),'g', tiempo_total, IC_maple_total,'r', corriente_Data(:, 1)  + 0.01467, corriente_Data(:, 2)/R3 + 6.6668099e-4,'b');
title('Corriente del capacitor (subamortiguado)')
xlabel('Tiempo (seg)')
ylabel('Corriente (A)')
legend('ODE','Maple','Experimental')
grid on;

% Inductor
corriente_inductor = valores_ode_total(:,2);
subplot(2, 2, 3);
plot(tiempo_total, corriente_inductor, 'g', tiempo_total, IL_maple_total, 'r', corriente2_Data(:, 1) + 0.0147, corriente2_Data(:, 2)/(R2-RL)  + 0.0012,'b');
title('Corriente del inductor (subamortiguado)')
xlabel('Tiempo (seg)')
ylabel('Corriente (A)')
legend('ODE','Maple','Experimental')
grid on;

voltaje_inductor = valores_ode_total(:,4) + (RL * valores_ode_total(:,2));
voltaje_inductor_maple = VL_maple_total + (RL * valores_ode_total(:,2));

subplot(2, 2, 4);
plot(tiempo_total, voltaje_inductor, 'g', tiempo_total, voltaje_inductor_maple, 'r', voltaje2_Data(:, 1) + 0.0147, voltaje2_Data(:, 2),'b');
title('Voltaje del inductor (subamortiguado)')
xlabel('Tiempo (seg)')
ylabel('Voltaje (V)')
legend('ODE','Maple','Experimental')
grid on;


%Diagrama de fases
figure(2);

% Gráfica de IL vs. VC
subplot(1,2,1);

plot(valores_ode_total(:,2), valores_ode_total(:,1), 'xg');
title('Diagrama de fases continua');
xlabel('IL');
ylabel('VC');
legend('Data ODE');

% Gráfica de IC vs. VL
subplot(1,2,2);
plot(valores_ode_total(:,3), voltaje_inductor, 'xg');
title('Diagrama de fases discontinua');
xlabel('IC');
ylabel('VL');
legend('Data ODE');
grid on;

% Errores

% VC
tiempo  = voltaje_Data(:, 1) + 0.01465;
voltaje = voltaje_Data(:, 2);
vcData  = [tiempo voltaje];
%V=-V;
condicion_inicial =  -R2*V/(R1+R2);
condicion_inicial_diff = 2*V/(C*(R1+R3));

tt = vcData(68:433,1);
VC_maple_ = condicion_inicial.*cosh((1/2)*tt*sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2)/(C*(R1+R3)*L)).*exp(-(1/2)*(C*R1*R2+C*R1*R3+C*R2*R3+L)*tt/(C*(R1+R3)*L))+(-cosh((1/2)*tt*sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2)/(C*(R1+R3)*L)).*exp(-(1/2)*(C*R1*R2+C*R1*R3+C*R2*R3+L)*tt/(C*(R1+R3)*L))+1)*R2*V/(R1+R2)+sinh((1/2)*tt*sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2)/(C*(R1+R3)*L)).*exp(-(1/2)*(C*R1*R2+C*R1*R3+C*R2*R3+L)*tt/(C*(R1+R3)*L))*(2*C*L*condicion_inicial_diff*R1*(R1+R2)+2*C*L*condicion_inicial_diff*R3*(R1+R2)+C*condicion_inicial*R1*R2*(R1+R2)+C*condicion_inicial*R1*R3*(R1+R2)+C*condicion_inicial*R2*R3*(R1+R2)-C*R1*R2^2*V-C*R2*R3*V*(R1+R2)+L*condicion_inicial*(R1+R2)-L*R2*V)/(sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2)*(R1+R2));
longitud = length(tt);
errorVmape = sum(abs(((vcData(68:433,2)+5) - (VC_maple_ + 5)) ./ (vcData(68:433,2)+5))) * 100 / longitud;
%figure(3)
%plot(tt, VC_maple_, tt, vcData(68:433,2));

% IC
tiempo = corriente_Data(:, 1) + 0.01467;
corriente = corriente_Data(:, 2)/R3 + 6.6668099e-4;
icData  = [tiempo corriente];
condicion_inicial = 2*V/(R1+R3);
condicion_inicial_diff = -2*V*(C*R1*R3+L)/(C*(R1+R3)^2*L);

tt = icData(68:433,1);
IC_maple_ = exp(-(1/2)*(C*R1*R2+C*R1*R3+C*R2*R3+L)*tt/(C*(R1+R3)*L)).*(condicion_inicial.*cosh((1/2)*tt*sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2)/(C*(R1+R3)*L))+sinh((1/2)*tt*sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2)/(C*(R1+R3)*L))*(2*condicion_inicial_diff*C*L*R1+2*condicion_inicial_diff*C*L*R3+C*condicion_inicial*R1*R2+C*condicion_inicial*R1*R3+C*condicion_inicial*R2*R3+L*condicion_inicial)/sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2));
longitud = length(tt);
errorImape = sum(abs(((icData(68:433,2)+5) - (IC_maple_ + 5)) ./ (icData(68:433,2)+5))) * 100 / longitud;
%figure(3)
%plot(tt, IC_maple_, tt,icData(68:433,2));

% IL
tiempo = corriente2_Data(:, 1) + 0.0147;
corriente = corriente2_Data(:, 2)/(R2-RL)  + 0.0012;
ilData  = [tiempo corriente];
condicion_inicial = -V/(R1+R2);
condicion_inicial_diff = 2*V*R3/((R1+R3)*L);

tt = ilData(270:733,1);
IL_maple_ = condicion_inicial.*cosh((1/2)*tt*sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2)/(C*(R1+R3)*L)).*exp(-(1/2)*(C*R1*R2+C*R1*R3+C*R2*R3+L)*tt/(C*(R1+R3)*L))+(-cosh((1/2)*tt*sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2)/(C*(R1+R3)*L)).*exp(-(1/2)*(C*R1*R2+C*R1*R3+C*R2*R3+L)*tt/(C*(R1+R3)*L))+1)*V/(R1+R2)+sinh((1/2)*tt*sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2)/(C*(R1+R3)*L)).*exp(-(1/2)*(C*R1*R2+C*R1*R3+C*R2*R3+L)*tt/(C*(R1+R3)*L))*(2*C*L*condicion_inicial_diff*R1*(R1+R2)+2*C*L*condicion_inicial_diff*R3*(R1+R2)+C*condicion_inicial*R1*R2*(R1+R2)+C*condicion_inicial*R1*R3*(R1+R2)+C*condicion_inicial*R2*R3*(R1+R2)-C*R1*R2*V-C*R3*V*(R1+R2)+L*condicion_inicial*(R1+R2)-L*V)/(sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2)*(R1+R2));
longitud = length(tt);
errorILmape = sum(abs(((ilData(270:733,2)+5) - (IL_maple_ + 5)) ./ (ilData(270:733,2)+5))) * 100 / longitud;
%figure(3)
%plot(tt, IL_maple_, tt, ilData(270:733,2));

% VL
tiempo = voltaje2_Data(:, 1) + 0.0147;
voltaje = voltaje2_Data(:, 2);
vlData  = [tiempo voltaje];
condicion_inicial = 2*V*R3/(R1+R3);
condicion_inicial_diff = -2*V*(C*R1*R2*R3+C*R1*R3^2+C*R2*R3^2-L*R1)/(C*(R1+R3)^2*L);

tt = vlData(270:733,1);
VL_maple_ = exp(-(1/2)*(C*R1*R2+C*R1*R3+C*R2*R3+L)*tt/(C*(R1+R3)*L)).*(condicion_inicial.*cosh((1/2)*tt*sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2)/(C*(R1+R3)*L))+sinh((1/2)*tt*sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2)/(C*(R1+R3)*L))*(2*condicion_inicial_diff*C*L*R1+2*condicion_inicial_diff*C*L*R3+C*condicion_inicial*R1*R2+C*condiciones_iniciales(4)*R1*R3+C*condicion_inicial*R2*R3+L*condicion_inicial)/sqrt(C^2*R1^2*R2^2+2*C^2*R1^2*R2*R3+C^2*R1^2*R3^2+2*C^2*R1*R2^2*R3+2*C^2*R1*R2*R3^2+C^2*R2^2*R3^2-4*C*L*R1^2-2*C*L*R1*R2-2*C*L*R1*R3-2*C*L*R2*R3+L^2));
VL_maple_total = VL_maple_ + (RL * ilData(270:733,2));
errorVLmape = sum(abs(((vlData(270:733,2)+5) - (VL_maple_total + 5)) ./ (vlData(270:733,2)+5))) * 100 / longitud;
%figure(3)
%plot(tt, VL_maple_total, tt, vlData(270:733,2));

format long;
errores = [errorVmape errorImape errorILmape errorVLmape];
disp('ERRORES VC, IC, IL, VL:')
disp(errores)
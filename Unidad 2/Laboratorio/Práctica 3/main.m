%Práctica 3 - Castro José, Hurtado Carlos.
clear all
close all
global Req R2 R3 Rm V C FREQ PERIODO;
%Parámetros del sistema

FREQ = 378.8; 
PERIODO = 1/FREQ;
Req = 1610;
R2 = 1020;
R3 = 2030;
Rm = 1030;
C = 0.1e-6;
V = 5.55;

%Variables de condiciones iniciales
vcIni = R2*V/(-R2-R3-Req); 
icIni = 2*R2*V/(R2*R3+R2*Req+R2*Rm+R3*Rm+Req*Rm);

%Variables de tiempo
tInicial=0;
tFinal=0;
tRango=[];
tTotal=[0];

condIni = [vcIni icIni];
%Vector para acumular los valores de x en el ODE
infODE = condIni;

%Variables con las condiciones iniciales de MAPLE
vcMaple = 0;
icMaple = 0;
%Vectores para acumular las C.I. de MAPLE
VcODE = [vcIni];
IcODE = [icIni];

 for i=1:2

    tRango = [tInicial tTotal(end)+PERIODO/2];     
    %condIni = [infODE(end,1) icIni];
    
    [t,x] = ode23('primerOrden',tRango,condIni);
    infODE = [infODE;x];
    tTotal = [tTotal;t];

    %Actualizar info MAPLE
    auxV = infODE(:,1);
    %disp(auxV);
    
    vcMaple = (((R2+R3+Req)*condIni(:,1)-R2*V)*exp(-(R2+R3+Req)*(t-tInicial)/(((R3+Req+Rm)*R2+Rm*(R3+Req))*C))+R2*V)/(R2+R3+Req);
    VcODE = [VcODE;vcMaple]; 
    icMaple = exp(-(R2+R3+Req)*(t-tInicial)/((R2*R3+R2*Req+R2*Rm+R3*Rm+Req*Rm)*C))*condIni(:,2);
    IcODE = [IcODE;icMaple];
    
    V = -V;
    vcIni = R2*V/(-R2-R3-Req);
    icIni = 2*R2*V/(R2*R3+R2*Req+R2*Rm+R3*Rm+Req*Rm);
    condIni = [vcIni icIni];
    tInicial = tTotal(end);
    tRango = [tTotal(end) tTotal(end)+PERIODO/2];  
       
    [t,x] = ode23('primerOrden',tRango,condIni);
    infODE = [infODE;x];
    tTotal = [tTotal;t];  
    
    vcMaple = (((R2+R3+Req)*condIni(:,1)-R2*V)*exp(-(R2+R3+Req)*(t-tInicial)/(((R3+Req+Rm)*R2+Rm*(R3+Req))*C))+R2*V)/(R2+R3+Req);
    VcODE = [VcODE;vcMaple]; 
    icMaple = exp(-(R2+R3+Req)*(t-tInicial)/((R2*R3+R2*Req+R2*Rm+R3*Rm+Req*Rm)*C))*condIni(:,2);
    IcODE = [IcODE;icMaple];
    
    V = -V;
    vcIni = R2*V/(-R2-R3-Req);
    icIni = 2*R2*V/(R2*R3+R2*Req+R2*Rm+R3*Rm+Req*Rm);
    
    condIni = [vcIni icIni];
    tInicial = tTotal(end);
 end

%Extracción de la data
dvoltaje = importdata('Data/Voltaje.CSV',',',18);
VoltajeData =  dvoltaje.data(:, 4:5);
dcorriente = importdata('Data/Corriente.CSV',',',18);
CorrienteData = dcorriente.data(:, 4:5);

%Gráficas Vc e Ic
figure(1);
plot((tTotal),infODE(:,1),'g',(tTotal),VcODE,'r', (VoltajeData(:,1))+2.8594e-3, VoltajeData(:,2)-2.32,'b');
title('Voltaje del capacitor');
xlabel('Tiempo (seg)');
ylabel('Voltaje (Vols)');
legend('ODE', 'Maple', 'Experimental');
grid on;

figure(2);
plot((tTotal),infODE(:,2),'g',(tTotal),IcODE,'r',(CorrienteData(:,1))+2.8594e-3,((CorrienteData(:,2))/Rm)+9.677e-7);
title('Corriente del capacitor');
xlabel('Tiempo (seg)');
ylabel('Voltaje (Vols)');
legend('ODE', 'Maple', 'Experimental');
grid on;

%Cálculo del error 
longitud = length(tTotal);
%errorVmae = mean(abs((VoltajeData(1:longitud,2) - VcODE)/ longitud))
%errorImae = mean(abs(((CorrienteData(1:longitud,1)/Rm) - IcODE) / longitud))
errorVmape = mean(abs((VoltajeData(1:longitud,2) - VcODE)./ VoltajeData(1:longitud,2)))/ longitud
errorImape = mean(abs(((CorrienteData(1:longitud,2)/Rm) - IcODE)./(CorrienteData(1:longitud,1)/Rm)))/ longitud



% t representa la variable indpendiente, 
% x representa el vector de variables de estados, 
% diff representa la primera derivada de las variables de estado

function dx = primerOrden(~, x) 
    %@brief Declaración de las variables Globales
    global Req R2 R3 Rm V C
 
    % Definición del modelo
    dx(1) = ((V-x(1))*R2-x(1)*R3-x(1)*Req)/((R2*R3+R2*Req+R2*Rm+R3*Rm+Req*Rm)*C);
    dx(2)= ((V-x(2)/C)*R2-x(2)*R3/C-x(2)*Req/C)/(R2*R3+R2*Req+R2*Rm+R3*Rm+Req*Rm);
    dx = dx';
end
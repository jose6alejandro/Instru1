function dx = ecuaciones(t,x)

global V R1 R2 R3 C L 

dx(1) = (V-x(1)-R1*x(2))/((R1+R3)*C);                                   % VC 
dx(2) = (-x(2)*R1*R2-x(2)*R1*R3-x(2)*R2*R3+x(1)*R1+V*R3)/((R1+R3)*L);   % IL 

dx(3) = (-x(3)/C-R1*x(4)/L)/(R1+R3);                                    % IC
dx(4) = (-x(4)*R1*R2/L-x(4)*R1*R3/L-x(4)*R2*R3/L+x(3)*R1/C)/(R1+R3);    % VL

dx=dx';
% Parámetros del sistema
m = 1;  % relación m1/m2
L1 = 1;  % longitud del primer péndulo
L2 = 2;  % longitud del segundo péndulo
x0 = [pi/3, 0, -pi/6, 0];  % condiciones iniciales [theta1, omega1, theta2, omega2]
tspan = [0, 10];  % intervalo de tiempo para la simulación

% Definición de la función del sistema de ecuaciones diferenciales
fg = @(t, x) [
    x(2);
    (cos(x(1) - x(3)) * (9.8 * sin(x(3)) / L1 - x(4)^2 * sin(x(1) - x(3))) - ...
    (L2 * ((m + 1) * 9.8 * sin(x(1)) / L2 + x(2)^2 * sin(x(1) - x(3))) / L1)) / ...
    (m + sin(x(1) - x(3))^2);
    x(4);
    (cos(x(1) - x(3)) * ((m + 1) * 9.8 * sin(x(1)) / L2 + x(2)^2 * sin(x(1) - x(3))) - ...
    ((m + 1) * L1 * (9.8 * sin(x(3)) / L1 - x(4)^2 * sin(x(1) - x(3))) / L2)) / ...
    (m + sin(x(1) - x(3))^2)
];


[t, x] = ode45(fg, tspan, x0);
% Gráfica de los resultados
figure;
subplot(2,1,1);
hold on
plot(t, x(:,1), 'r'); % Graficar theta1 en rojo
plot(t, x(:,3), 'b'); % Graficar theta2 en azul
hold off
grid on
legend('\theta_1', '\theta_2');
xlabel('Tiempo (s)');
ylabel('Ángulo (rad)');
title('Ángulos del péndulo doble');

subplot(2,1,2);
hold on
plot(t, x(:,2), 'r'); % Graficar omega1 en rojo
plot(t, x(:,4), 'b'); % Graficar omega2 en azul
hold off
grid on
legend('\omega_1', '\omega_2');
xlabel('Tiempo (s)');
ylabel('Velocidad angular (rad/s)');
title('Velocidades angulares del péndulo doble');

% Cálculo de la energía inicial 
Einicial = (m + 1) * L1 * x0(2)^2 / (2 * L2) + L2 * x0(4)^2 / (2 * L1) + ...
    (x0(2) * x0(4)) * cos(x0(1) - x0(3)) - (m + 1) * 9.8 * cos(x0(1)) / L2 - 9.8 * cos(x0(3)) / L1;

E = (m + 1) * L1 * x(:,2).^2 / (2 * L2) + L2 * x(:,4).^2 / (2 * L1) + ...
    (x(:,2) .* x(:,4)) .* cos(x(:,1) - x(:,3)) - (m + 1) * 9.8 * cos(x(:,1)) / L2 - 9.8 * cos(x(:,3)) / L1;

% Calcular el error máximo relativo de energía
error_maximo = max(abs(100 * (Einicial - E) / Einicial));
disp(['Error máximo relativo de energía: ', num2str(error_maximo), '%']);

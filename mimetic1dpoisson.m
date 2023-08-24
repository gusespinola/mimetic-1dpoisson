% Comentar cuando modo funcion = off
%function v = test20230331(h)
% Sintaxis
% v = test20230331(h);
% Parametros internos
% N              cantidad de nodos                       
% a0             nodo inicial
% b0             nodo final
% h = (b0-a0)/N  tamanho de celda 1D
% lambda         parametro de la ecuacion dif.
% alpha          parametro de la cond. de Robin
% beta           parametro de la cond. de Robin
% epsilon        vertical shift para mejor visualizacion

a0 = 0;
b0 = 1;

% % Test 1 (comentar si Test 2 = ON)
% N = 100; h = (b0-a0)/N;

% Test 2
h = 0.005; % cambiar a 0.20, 0.10, 0.05
N = (b0 - a0) / h;

% Condiciones de Robin
lambda = -1;
alpha = -exp(lambda);
beta = (exp(lambda) - 1)/ lambda;

% Funciones de prueba (ver mimeticos_tesis_2022.pdf, cap 3)
% Funcion RHS
F  = @(x) lambda*lambda*exp(lambda*x)/(exp(lambda)-1);
% Funcion solucion de la ED
f  = @(x) (exp(lambda*x)-1)/(exp(lambda)-1);

% Matriz gradiente G
G = zeros(N+1, N+2);
G(1,1:3) = [-8/3 3 -1/3];
for i=2:N
    G(i, i:i+1) = [-1 1];
end
G(N+1, N:N+2) = [1/3 -3 8/3];
G = (1/h)*G;

% Matriz divergencia D
D = zeros(N, N+1);
for i = 1:N
    D(i, i:i+1) = [-1 1];
end
D = (1/h)*D;

% Matriz L
L = zeros(N+2, N+2);
L(2:N+1, :) = D*G;

% Matriz divergencia ampliada Dt, o "D tilde";
Dt = zeros(N+2, N+1);
Dt(2:N+1, :) = D;

% Matriz del producto escalar P1 (trapecio)
P1 = zeros(N+1);
P1(1,1) = 1/2;
P1(N+1, N+1) = 1/2;
P1(2:N, 2:N) = eye(N-1);

% Matriz del producto escalar P2 (simpson?)
P2 = zeros(N+1);
P2(1,1) = 3/8;
P2(2,2) = 9/8;
P2(N, N) = 9/8;
P2(N+1, N+1) = 3/8;
P2(3:N-1, 3:N-1) = eye(N-3);

% Matriz B1
B1 = zeros(N+2, N+1);
B1(1,1)= -1;
B1(N+2, N+1) = 1;
B11 = B1*G;

% Matriz B2
B2 = Dt + G' * P2;
B2 = h*B2;
B22 = B2*G;

% RHS b
b = zeros(N+2,1);
b(1,1) = -1;
b(N+2,1) = 0;
xi = i*h;
for i = 2:N+1
    b(i,1) = F((i-1-1/2)*h);
end

% Matriz A
A = zeros(N+2, N+2);
A(1,1) = 1;
A(N+2, N+2) = 1;

% Matriz LHS # 1
A1 = alpha*A + beta*B1*G + L;
cond=condest(A1);

% Matriz LHS # 2
A2 = alpha*A + beta*B2*G + L;

% Sistema lineal # 1
%        [L,U] = ilu(sparse(A1),struct('type','ilutp','droptol',1e-6));
%        M=L*U; 
%        A1=M\A1; b=M\b;
%x1 = A1\b;
%x1 = inv(A1)*b;
%m=N/2;

 %tStart=tic;
%[x1,flag1,resgmres1,ITER1,resvec1] = gmres(A1,b,10,1e-9,1000);
[x1,flag1,resgmres1,ITER1,resvec1] = bicgstab(A1,b,1e-9,1000);
%time =toc(tStart)

% Sistema lineal # 2
%x2 = A2\b;
%x2 = inv(A2)*b;


% Shift x2 hacia arriba para mejor visualizacion
% for i = 1:size(x2,1)
%     x2(i,1) = x2(i,1); %+ epsilon; % quitar el epsilon si es necesario
% end

% Grid: xcb = valores de x en los nodos y en la frontera
xcb = zeros(1, N+2);
xcb(1,1) = a0;
xcb(1,N+2) = b0;
for i = 2:N+1
    xcb(1,i) = a0 + (i-1-1/2)*h;
end
xcb = xcb';

% Plot
figure(1)
ylim([0 1])
plot(xcb, f(xcb), xcb, x1)
legend('f(xcb)', 'x1')

% Medir error RMSE, ver Tesis
% e1 = xcb(2:N+1) - x1(2:N+1);
% e2 = xcb(2:N+1) - x2(2:N+1);
%e1 = xcb - x1;
f_xcb=f(xcb);
e1 = f_xcb - x1;
error=norm(e1);
error_relativo=norm(e1)/norm(f_xcb);
%e2 = xcb - x2;
%e11 = e1.*e1;
%e22 = e2.*e2;
%av1 = mean(e11);
%av2 = mean(e22);
%rmse1 = sqrt(e11);
ecm = immse(f_xcb,x1); %Error cuadratico medio
%rmse2 = sqrt(av2);

% Imprimir resultados
% fprintf('m = %f\n', m)

fprintf('h = %f\n', h)
%fprintf('av1 = %f\n', av1)
%fprintf('av2 = %f\n', av2)
%fprintf('rmse1 = %f\n', rmse1)
fprintf('ecm = %f\n', ecm)
fprintf('error = %f\n', error)
fprintf('error_relativo = %f\n', error_relativo)


fprintf('cond = %f\n', cond)

%v = [h av1 av2 rmse1 rmse2];

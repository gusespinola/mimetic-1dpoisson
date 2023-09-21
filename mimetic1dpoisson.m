clear all
clc
%% Problem parameters - 1D Poisson Equation                       
% a0             initial value of input x
% b0             final value of input x
% h              step size
% alpha, beta    Robin constraint parameters

a0 = 0;
b0 = 1;
h = 0.1; % Change to: 0.20, 0.10, 0.05
N = (b0-a0)/h; % Number of nodes

% Robin boundary conditions
lambda = -1;
alpha = -exp(lambda);
beta = (exp(lambda) - 1)/ lambda;

% Functions for 1D Poisson Equation
% \nabla^2 u(x) = F(x)
% see \cite{hernandez2007large} for more information
% RHS function of the 1D Poisson Equation
F  = @(x) lambda*lambda*exp(lambda*x)/(exp(lambda)-1);
% Analytical solution of the 1D Poisson Equation
f  = @(x) (exp(lambda*x)-1)/(exp(lambda)-1);

% Gradient opeator G
G = zeros(N+1, N+2);
G(1,1:3) = [-8/3 3 -1/3];
for i=2:N
    G(i, i:i+1) = [-1 1];
end
G(N+1, N:N+2) = [1/3 -3 8/3];
G = (1/h)*G;

% Divergence operator D
D = zeros(N, N+1);
for i = 1:N
    D(i, i:i+1) = [-1 1];
end
D = (1/h)*D;

% Laplace operator L
L = zeros(N+2, N+2);
L(2:N+1, :) = D*G;

% Expanded Divergence operator Dt
Dt = zeros(N+2, N+1);
Dt(2:N+1, :) = D;

% Weighted inner product matrix P1
P1 = zeros(N+1);
P1(1,1) = 1/2;
P1(N+1, N+1) = 1/2;
P1(2:N, 2:N) = eye(N-1);

% Weighted inner product matrix P2
P2 = zeros(N+1);
P2(1,1) = 3/8;
P2(2,2) = 9/8;
P2(N, N) = 9/8;
P2(N+1, N+1) = 3/8;
P2(3:N-1, 3:N-1) = eye(N-3);

% Boundary operator B1
B1 = zeros(N+2, N+1);
B1(1,1)= -1;
B1(N+2, N+1) = 1;
B11 = B1*G;

% Boundary operator B2
B2 = Dt + G' * P2;
B2 = h*B2;
B22 = B2*G;

% RHS vector 'b' for the linear system Ai*x = b, i=1,2
b = zeros(N+2,1);
b(1,1) = -1;
b(N+2,1) = 0;
xi = i*h;
for i = 2:N+1
    b(i,1) = F((i-1-1/2)*h);
end

% Matrix A
A = zeros(N+2, N+2);
A(1,1) = 1;
A(N+2, N+2) = 1;

% Matriz LHS # 1 by Support Operator Method
A1 = alpha*A + beta*B1*G + L;
cond=condest(A1);

% Matriz LHS # 2 by Castillo-Grone Method
A2 = alpha*A + beta*B2*G + L;
condA2=condest(A2);

% Linear system # 1
%        [L,U] = ilu(sparse(A1),struct('type','ilutp','droptol',1e-6));
%        M=L*U; 
%        A1=M\A1; b=M\b;
%x1 = A1\b;
%x1 = inv(A1)*b;
%m=N/2;

%% Preconditioners
%        [L,U] = ilu(A);
        setup.type = 'ilutp';
         setup.droptol = 1e-6;
         setup.udiag=1;
       [L,U] = ilu(sparse(A1),setup);
         M=L*U; 
%         clear L U;
         A1=M\A1; b=M\b;


% Jacobi
% CondA=condest(A)
% tStart=tic;
% M_jacobi=(diag(diag(A1)))^-1;
% A1=M_jacobi*A1;
% b=M_jacobi*b;
% cond=condest(A1);

% SOR
% CondA=condest(A1)
% n=size(A);
% tStart=tic;
% omega=1;
% D=diag(diag(A1));
% L=-(tril(A1)-D);
% M_sor=(1/omega)*D-L;
% A1=M_sor\A1;
% b=M_sor\b;

%% Iterative methods
 %tStart=tic;
%[x1,flag1,resgmres1,ITER1,resvec1] = gmres(A1,b,10,1e-6,1000);
%x1 = gmres(A1,b,10,1e-9,1000);
[x1,flag1,resgmres1,ITER1,resvec1] = bicgstab(A1,b,1e-6,1000);
%time =toc(tStart)

%% Adaptive Iterative methods - UNDER CONSTRUCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Par�metros globales
% itermax = 1000; %2000 funciona
% tol = 1e-06;
% end_count = 1;
% % Par�metro de Adaptive-LGMRESE     
% eps0 = 0.5;
% alpha = 3;
% Name_Matrix_1 = 'A_1 u=b';
% Name_Matrix_2 = 'A_2 u=b';
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % LGMRES-E(m,d,l)
% m0 = 30;
% d0 = 2;
% l0 = 2;
% color = 'r-';
% print = 0;
% temp = zeros(end_count, 1);
% for count = 1:end_count
%     if count == end_count
%         print = 1;
%     end
%     [vec5, vny, xx] = adaptive_LGMRESE(A1, b, m0, d0, l0, itermax, tol, alpha, eps0, color, print, Name_Matrix_1);
%     temp(count)=vec5(1);
% end
% sum_s= vec5(3);
% t_prom5 = mean(temp);
% t_std_dev5 = std(temp);
% restart = vec5(2);
% metrics_LGMRESE = [t_prom5 t_std_dev5 vec5(2) vec5(3)];

%% Parameters
% alpha=-3; %2
% delta=5; %0.8
% opts_tol=1e-6;
% itermax=1000;
% p = 1;
    
% %%       %PD-GMRES(m)
% % The original idea was to compute the average execution time,
% % we may discuss if this is still necessary 
% for i=1:p
%     color_pd_gmres='b';
%     mPD=30;
%     %alpha=2;
%     %rootFolder = fileparts(pwd); % go to the root folder
%     %srcFolder = fullfile(rootFolder, 'src'); % enter the data folder
%     %cd(srcFolder)
%     [time, logres_pd_gmres, u3]=pd_gmres(A1,b, mPD, alpha, delta,itermax, Name_Matrix_1);
%     %sol7(size(sol7,1)+1,:)= [time];
% end

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
figure
p = plot(xcb, f(xcb), 'r-x', ...
    xcb, x1, 'b--+');
p(1).LineWidth = 2;
p(2).LineWidth = 2;
xlim([0 1])
ylim([0 1])
legend('Analytical solution', 'Numerical solution', 'Location', 'southeast')
title('Numerical solution v. analytical solution',...
    'FontWeight','bold')
xlabel('x_i in [0,1]')
ylabel('u_i = u(x_i)')

% Medir error RMSE, ver Tesis
%e1 = xcb(2:N+1) - x1(2:N+1);
% e2 = xcb(2:N+1) - x2(2:N+1);
%e1 = xcb - x1;
f_xcb=f(xcb);
e1 = f_xcb - x1;
%error=norm(e1);
%error_relativo=norm(e1)/norm(f_xcb);
%e2 = xcb - x2;
e11 = e1.*e1;
%e22 = e2.*e2;
av1 = mean(e11);
%av2 = mean(e22);
ecm = sqrt(av1);


%ecm = immse(f_xcb,x1); %Error cuadrático medio


%rmse2 = sqrt(av2);

% Imprimir resultados
% fprintf('m = %f\n', m)

fprintf('h = %f\n', h)
%fprintf('av1 = %f\n', av1)
%fprintf('av2 = %f\n', av2)
%fprintf('rmse1 = %f\n', rmse1)
fprintf('ecm = %f\n', ecm)
%fprintf('error = %f\n', error)
%fprintf('error_relativo = %f\n', error_relativo)


fprintf('cond = %f\n', cond)

%v = [h av1 av2 rmse1 rmse2];
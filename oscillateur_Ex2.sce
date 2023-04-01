clear;
k=1; m=1;
A = [0     1;
    (-k/m) 0];

function xdot = f(t, x)
    xdot = A*x;
endfunction

pos0 = 0; 
u0 = 1; 
x0 = [pos0; u0]; // Donnees initiale (colonne)
t0 = 0;
N = 1000; 
T=10*%pi; 
tau = (T-t0)/N;
td = t0 : tau : T;

function Xsol = RK2(x0, t0, T, N, f)
    x_prev=x0;
    x_curr=[0;0];
    x_curr_chap=[0;0];
    Xsol = [];
    for i = 1:N+1
        x_curr_chap=x_prev + tau*f(i*tau, x_prev);
        x_curr = x_prev + (tau/2)*(f(i*tau, x_prev)+f((i+1)*tau, x_curr_chap));
        x_prev = x_curr;
        Xsol(:,i) = x_curr;
    end;
endfunction


Xsol=RK2(x0,t0,T,N,f);

figure(1); clf;
//plot(td, Xsol(1,:), '.-');ylabel('x'); xgrid;
subplot(2,1,1), plot(td, Xsol(1,:), '.-');ylabel('x'); xgrid;
subplot(2,1,2), plot(td, Xsol(2,:), '.-');ylabel('u'); xgrid;
//pause
figure(2); clf;
plot(Xsol(1,:), Xsol(2,:), '-'); mtlb_axis('equal');
comet(Xsol(1,:), Xsol(2,:), 0.01);
//pause
figure(3); clf;
comet(Xsol(1,:), 0*Xsol(1,:), 0.01);
// Calcul de l'energie total calculée en cours du temps
//
xt = Xsol(1,:);
ut = Xsol(2,:);
Et = 0.5*k * xt.^2 +0.5*m * ut.^2;
figure(4); clf;
mtlb_axis([t0 T 0 1]);
plot(td, Et);

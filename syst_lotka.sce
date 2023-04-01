// Systeme de Lotka-Volterra
//
clear;

function xdot = Lotka(t, xvec)
    x = xvec(1);
    y = xvec(2);
    xdot(1) = (1-y)*x;
    xdot(2) = (x-1)*y;
endfunction

x0 = [3;1]; // vecteur colonne
t0 = 0;
T = 80;
N = 6000;
tau = T / N;
td = t0 : tau : T;
Xsol = ode(x0, t0, td, Lotka);
//
figure(1); clf;
//plot(td, Xsol(1,:), '.-');ylabel('x'); xgrid;
subplot(2,1,1), plot(td, Xsol(1,:), '.-');ylabel('x'); xgrid;
subplot(2,1,2), plot(td, Xsol(2,:), '.-');ylabel('y'); xgrid;
//pause
figure(2); clf;
plot(Xsol(1,:), Xsol(2,:), '-'); mtlb_axis('equal');
comet(Xsol(1,:), Xsol(2,:), 0.01); xgrid;
//
xt = Xsol(1,:);
yt = Xsol(2,:);
eta = xt - log(xt) + yt - log(yt);
figure(3); clf;
mtlb_axis([t0 T 2 3]);
plot(td, eta, '.-');

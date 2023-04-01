clear;
k=1; m=1;
A = [0     1;
    (-k/m) 0];

function xdot = LotkaVolterra(t, xvec)
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


Xsol=RK2(x0,t0,T,N,LotkaVolterra);

Xsol2 = ode(x0, t0, td, LotkaVolterra);

figure(1); clf;
//plot(td, Xsol(1,:), '.-');ylabel('x'); xgrid;
subplot(2,1,1), plot(td, Xsol(1,:), 'r.-' );ylabel('x'); xgrid;
subplot(2,1,2), plot(td, Xsol(2,:), 'r.-');ylabel('y'); xgrid;

subplot(2,1,1), plot(td, Xsol2(1,:), 'g.->');ylabel('x'); xgrid;
subplot(2,1,2), plot(td, Xsol2(2,:), 'g.->');ylabel('y'); xgrid;

//pause
figure(2); clf;
plot(Xsol(1,:), Xsol(2,:), 'r-'); mtlb_axis('equal');
comet(Xsol(1,:), Xsol(2,:), 0.01); xgrid;

plot(Xsol2(1,:), Xsol2(2,:), 'g-'); mtlb_axis('equal');
comet(Xsol2(1,:), Xsol2(2,:), 0.01); xgrid;
//
xt = Xsol(1,:);
yt = Xsol(2,:);
eta = xt - log(xt) + yt - log(yt);

xt2 = Xsol2(1,:);
yt2 = Xsol2(2,:);
eta2 = xt2 - log(xt2) + yt2 - log(yt2);

figure(3); clf;
mtlb_axis([t0 T 2 3]);
plot(td, eta, 'r.-');

mtlb_axis([t0 T 2 3]);
plot(td, eta2, 'g.-');

clear;

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

//Xsol = ode(x0,t0,td,LotkaVolterra);

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

Xsol = RK2(x0,t0,T,N,LotkaVolterra);


X = Xsol(:, 1:$-1);
Y = Xsol(:, 2:$);

Arond = Y * X' * inv(X * X');
[P, diagevals] = spec(Arond);

lamdak = zeros(size(x0)(1),size(x0)(1))

for i=1:size(x0)(1)
    lamdak(i,i) = log(diagevals(i,i))/tau;
end

x_curr = [0 0];
x_prev = x0;
Xsol_appr = []
for i = 1:N+1
    x_curr = Arond*x_prev;
    Xsol_appr(:,i) = x_curr;
    x_prev = x_curr;
end;







figure(1); clf;
//plot(td, Xsol(1,:), '.-');ylabel('x'); xgrid;
subplot(2,1,1), plot(td, Xsol(1,:), 'r.-');ylabel('x'); xgrid;
subplot(2,1,2), plot(td, Xsol(2,:), 'r.-');ylabel('y'); xgrid;

subplot(2,1,1), plot(td, Xsol_appr(1,:), 'g.-');ylabel('x'); xgrid;
subplot(2,1,2), plot(td, Xsol_appr(2,:), 'g.-');ylabel('y'); xgrid;
//pause
figure(2); clf;
plot(Xsol(1,:), Xsol(2,:), 'r-'); mtlb_axis('equal');
comet(Xsol(1,:), Xsol(2,:), 0.01); xgrid;

plot(Xsol(1,:), Xsol_appr(2,:), 'g-'); mtlb_axis('equal');
comet(Xsol(1,:), Xsol_appr(2,:), 0.01); xgrid;












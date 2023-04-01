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

x_point_chap = [];
for i = 3:N+1
    x_point_chap(:,i) = (3/(2*tau))*Xsol(:,i) ...
                        - (2/(tau))*Xsol(:,i-1) ...
                        + (1/(2*tau))*Xsol(:,i-2);
end;

phix = [Xsol(1,:); Xsol(2,:); Xsol(1,:).*Xsol(2,:)];

X = phix(:,3:$);
Y = x_point_chap(:,3:$);

A = Y * X' * inv(X * X');

function xdot = f(t, xvect)
    x = xvect(1);
    y = xvect(2);
    phix = [x;y;x*y];
    xdot = A*phix;
endfunction

Xsol_appr =  ode(x0, t0, td, f);

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

plot(Xsol_appr(1,:), Xsol_appr(2,:), 'g-'); mtlb_axis('equal');
comet(Xsol_appr(1,:), Xsol_appr(2,:), 0.01); xgrid;












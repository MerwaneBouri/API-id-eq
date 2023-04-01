
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

X = Xsol(:,1:$-1);
Y = Xsol(:,2:$);

B = Y * X' * inv(X * X');
E = B - expm(B*tau);

[P, diagevals] = spec(B);
lamda1 = log(diagevals(1,1))/tau;
lamda2 = log(diagevals(2,2))/tau;
//lamdak = log(diagevals)/tau

lamdak = zeros(2,2)
lamdak(1,1) = lamda1;
lamdak(2,2) = lamda2;

C = P * lamdak * inv(P);


x_curr = [0 0];
x_prev = x0;
Xsol_appr = []
for i = 1:N+1
    x_curr = B*x_prev;
    Xsol_appr(:,i) = x_curr;
    x_prev = x_curr;
end;



xt = Xsol_appr(1,:);
ut = Xsol_appr(2,:);
Et = 0.5*k * xt.^2 +0.5*m * ut.^2;
figure(4); clf;
mtlb_axis([t0 T 0 1]);
plot(td, Et);



figure(2); clf;
//plot(td, Xsol(1,:), '.-');ylabel('x'); xgrid;
subplot(2,1,1), plot(td, Xsol_appr(1,:), '.-');ylabel('x'); xgrid;
subplot(2,1,2), plot(td, Xsol_appr(2,:), '.-');ylabel('u'); xgrid;




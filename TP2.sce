clear;
//système ocillant multidimmensionnel 
r=6

K = 2*eye(r,r) - diag(ones(r-1,1),-1) -diag(ones(r-1,1),+1)

A=zeros(2*r,2*r)

A(1:r,r+1:$)=eye(r,r)

A(r+1:$,1:r)=-K

[P, diagevals]=spec(A)

function xdot=oscildim(t,x)
    xdot = A*x
endfunction

a0=zeros(r,1)
b0=zeros(r,1)
b0(1)=1

x0=[a0;b0]

t0=0; T = 120; N=1000;
tau = (T-t0)/N;
td = t0 : tau : T;
Xsol = ode(x0, t0, td, oscildim);


X = Xsol(:, 1:$-1);
Y = Xsol(:, 2:$);

Arond = Y * X' * inv(X * X');

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
for i=1:r
    subplot(r,1,i), plot(td, Xsol_appr(i,:), 'r.-' );//ylabel('x %i'); xgrid;
    subplot(r,1,i), plot(td, Xsol(i,:), 'g.-' );ylabel('xi'); xgrid;
end
figure(3); clf;
plot(Xsol_appr(1,:), Xsol_appr(7,:), 'r-'); mtlb_axis('equal');
comet(Xsol_appr(1,:), Xsol_appr(7,:), 0.01); xgrid;


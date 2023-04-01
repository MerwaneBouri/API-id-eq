m1 = 1;
m2 = 0.5;
k1 = 1;
k2 = 1;
k3 = 0.5;
N = 1000;
t0 = 0;
T = 80;
A = [0  1    0   0;
    (-(k1+k2)/m1)   0   (k2/m1) 0 ;
    0   0   0   1;
    (k2/m2) 0   (-(k2+k3)/m2)   0]; 

C =[0 ; -k2/m1; 0; ((k2+k3)/m2)];

tau = (T-t0)/N;
td = t0 : tau : T;

function xdot = f(t, x)
    //xdot = A*x+C;
    xdot = A*x;
endfunction

x0=[0;1;0;0]

Xsol =  ode(x0, t0, td, f);

X = Xsol(:, 1:$-1);
Y = Xsol(:, 2:$);

Arond = Y * X' * inv(X * X');
E = Arond - expm(Arond*tau);

[P, diagevals] = spec(Arond);

lamdak = zeros(size(x0)(1),size(x0)(1))

for i=1:size(x0)(1)
    lamdak(i,i) = log(diagevals(i,i))/tau;
end


B = P * lamdak * inv(P);

x_curr = [0 0];
x_prev = x0;
Xsol_appr = []
for i = 1:N+1
    x_curr = Arond*x_prev;
    Xsol_appr(:,i) = x_curr;
    x_prev = x_curr;
end;

Em = 1/2*k1*Xsol_appr(1,:)^2 ...
    + 1/2*k3*Xsol_appr(3,:)^2 ...
    + 1/2*k2*(Xsol_appr(3,:)-Xsol_appr(1,:))^2 ...
    + 1/2*m1*(Xsol_appr(2,:))^2 ...
    + 1/2*m2*(Xsol_appr(4,:))^2;
   
/*    
figure(1); clf;
//plot(td, Xsol(1,:), '.-');ylabel('x'); xgrid;
subplot(2,1,1), plot(td, Em, 'r.-' );ylabel('x'); xgrid;
*/

figure(1); clf;
er=4
//plot(td, Xsol(1,:), '.-');ylabel('x'); xgrid;
subplot(2,1,1), plot(td, Xsol(1,:), 'r.-' );ylabel('x1'); xgrid;
subplot(2,1,2), plot(td, Xsol(2,:), 'r.-');ylabel('u1'); xgrid;

subplot(2,1,1), plot(td, Xsol_appr(1,:), 'g.->');ylabel('x'); xgrid;
subplot(2,1,2), plot(td, Xsol_appr(2,:), 'g.->');ylabel('y'); xgrid;


//pause
figure(2); clf;
plot(Xsol(1,:), Xsol(2,:), 'r-'); mtlb_axis('equal');
comet(Xsol(1,:), Xsol(2,:), 0.01); xgrid;

plot(Xsol_appr(1,:), Xsol_appr(2,:), 'g-'); mtlb_axis('equal');
comet(Xsol_appr(1,:), Xsol_appr(2,:), 0.01); xgrid;






figure(1); clf;
//plot(td, Xsol(1,:), '.-');ylabel('x'); xgrid;
subplot(2,1,1), plot(td, Xsol(1,:), '.-');ylabel('x1'); xgrid;
subplot(2,1,2), plot(td, Xsol(2,:), '.-');ylabel('u1'); xgrid;

figure(2); clf;
subplot(2,1,1), plot(td, Xsol(3,:), '.-');ylabel('x2'); xgrid;
subplot(2,1,2), plot(td, Xsol(4,:), '.-');ylabel('u2'); xgrid;

//vitesse
figure(3); clf;
plot(Xsol(1,:), Xsol(2,:), 'r-'); mtlb_axis('equal');
comet(Xsol(1,:), Xsol(2,:), 0.01); xgrid;

plot(Xsol(3,:), Xsol(4,:), 'r-'); mtlb_axis('equal');
comet(Xsol(3,:), Xsol(4,:), 0.01); xgrid;
/*
figure(3); clf;
comet(Xsol(1,:), 0*Xsol(1,:), 0.01);

comet(Xsol(3,:), 0*Xsol(3,:), 0.01);
*/



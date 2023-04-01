sigma = 10;
rho = 28;
beta_ = 8/3;

function xdot = Lorenz(t, xvec)
    x = xvec(1);
    y = xvec(2);
    z = xvec(3);
    xdot(1) = sigma*(y-x);
    xdot(2) = rho * x - y - x * z;
    xdot(3) = x * y - beta_ * z;
endfunction

T =75;
N = 15000;
tau = T/N;
t0 = 0;

x0 = [0; 1; 1.05];
td = t0 : tau : T;

Xsol = ode(x0,t0,td,Lorenz);
/*
figure(1); clf;
comet3d(Xsol(1,:),Xsol(2,:),Xsol(3,:));
*/
x_point_chap = [];
for i = 3:N+1
    x_point_chap(:,i) = (3/(2*tau))*Xsol(:,i) ...
                        - (2/(tau))*Xsol(:,i-1) ...
                        + (1/(2*tau))*Xsol(:,i-2);
end;

phix = [Xsol(1,:); Xsol(2,:); Xsol(3,:); Xsol(1,:).*Xsol(2,:); Xsol(1,:).*Xsol(3,:);Xsol(2,:).*Xsol(3,:)];

X = phix(:,3:$);
Y = x_point_chap(:,3:$);

A = Y * X' * inv(X * X');

function xdot = f(t, xvect)
    x = xvect(1);
    y = xvect(2);
    z = xvect(3);
    phix = [x;y;z;x*y;x*z;y*z];
    xdot = A*phix;
endfunction

Xsol_appr =  ode(x0, t0, td, f);

figure(1); clf;
comet3d(Xsol(1,:),Xsol(2,:),Xsol(3,:));

figure(2); clf;
comet3d(Xsol_appr(1,:),Xsol_appr(2,:),Xsol_appr(3,:));

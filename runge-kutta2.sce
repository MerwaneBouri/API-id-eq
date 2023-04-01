


function xdot = f(t, x)
    xdot = A*x;
endfunction

function X = RK2(x0, t0, T, N, f)
	tau = (T-t0)/N
	x_prev=x0
	x_curr=0
	x_curr_chap=0
	Xsol = []
	for i = 1:N
		x_curr_chap=x_prev + tau*f(i*tau, x_prev)
		x_curr = x_prev + (tau/2)*(f(i*tau, x_prev)+f((i+1)*tau, x_curr_chap))
		Xsol(:,i) = x_curr
endfunction

x0 = [pos0; u0]; // Donnees initiale (colonne)
t0 = 0;
N = 4000; T=90; 
Xsol=RK2(x0,t0,T,N,f)
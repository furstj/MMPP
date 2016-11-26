clear all;

% Reseni smisene ulohy pro rovnici u_t + (u^2/2)_x = 0 pomoci MKO

% Delka intervalu a pocet bunek
L = 1;
n = 400;

dx = L / n;
x = linspace(dx/2, L-dx/2, n); 

u(1:n)  = 0;
un(1:n) = 0;

% Pocatecni podminka
for i = 1:n
  if (x(i)>0.25 && x(i)<0.5)
     u(i) = 1-cos(2*pi*(x(i)/0.25-1));
  end
end

plot(x, u); axis([0 1 -0.5 2.5]);
disp("Stiskni enter pro pokracovani"); pause;

dt = 0.1*dx;

t = 0;

for iter = 1:n

    % Okrajova podminka na leve casti intervalu
    f(1) = 0;
    f(n+1) = u(n)^2/2;

    % Numericky tok f
    for i = 2:n
	a = (u(i-1)+u(i))/2;
	q = a^2*dt/dx;
	%d = abs(a)-q;
	d = 25*abs(u(i)-u(i-1)) / (abs(u(i))+abs(u(i-1)) + 1.e-6) * (abs(a)-q);
	f(i) = a/2*(u(i-1)+u(i)) - q/2*(u(i)-u(i-1));
	f(i) = f(i) - d/2*(u(i)-u(i-1));
    end

    un(1:n) = u(1:n) - dt/dx * (f(2:n+1) - f(1:n));
    t = t + dt;

    u = un;

    if (mod(iter,10)==0)
      plot(x, u); axis([0 1 -0.5 2.5]);
      pause(1);
    end
end

t

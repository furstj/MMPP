%clear all;

% Reseni smisene ulohy pro rovnici u_t + a*u_x = 0 pomoci MKO
a = 1;

% Delka intervalu a pocet bunek
L = 1;
n = 100;

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

dt = 0.8*dx;

for iter = 1:(n/2)

    % Okrajova podminka na leve casti intervalu
    f(1) = 0;
    f(2) = 0;
    f(n+1) = a*u(n-1);
    f(n+1) = a*u(n);

    % Numericky tok f
    for i = 3:n-1
	fl = a*u(i-1);
	fr = a*u(i);
	flw = (fr+fl)/2 - a/2*dt/dx*(fr-fl);
	fup = (fr+fl)/2 - abs(a)/2*(u(i)-u(i-1));

	if ( (u(i)-u(i-1))*(u(i+1)-u(i)) <= 0) 
	   phi = 0;
	else
	  r = (u(i-1)-u(i-2)) / (u(i)-u(i-1));
	  phi = min(1, 2*r);
	end
	f(i) = fup + phi*(flw-fup);
    end

    un(1:n) = u(1:n) - dt/dx * (f(2:n+1) - f(1:n));

    u = un;

    if (mod(iter,10)==0)
      plot(x, u); axis([0 1 -0.5 2.5]);
      pause(1);
    end
end




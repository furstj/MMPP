clear all;

% Reseni smisene ulohy pro rovnici u_t + a*u_x = 0 pomoci MKO
a = 1;

% Delka intervalu a pocet bunek
L = 1;

% Cas ukonceni vyoctu
T = 0.25;


for n = 100:100:1000
  
  disp(n)

  dx = L / n;
  x = dx/2:dx:(L-dx/2);

  u(1:n)  = 0;    % aktualni hodnota reseni
  un(1:n) = 0;    % nova hodnota reseni
  
  % Pocatecni podminka
  for i = 1:n
    if (x(i)>0.25 && x(i)<0.5)
      u(i) = 1-cos(2*pi*(x(i)/0.25-1));
    end
  end

  % Presne reseni
  uex(1:n) = 0;
  for i = 1:n
    if (x(i)>0.25+a*T && x(i)<0.5+a*T)
      uex(i) = 1-cos(2*pi*((x(i)-a*T)/0.25-1));
    end
  end

  dt = 0.8*dx;
  t = 0;

  while (t<T) 

    if (t+dt>=T)
      dt = T - t;
      t = T;
    else
      t = t + dt;
    end

    % Okrajova podminka na leve casti intervalu
    f(1) = 0;
    f(n+1) = a*u(n);

    % Numericky tok f
    for i = 2:n
	f(i) = a/2*(u(i-1)+u(i)) - a^2/2*dt/dx*(u(i)-u(i-1));
    end
    
    un(1:n) = u(1:n) - dt/dx * (f(2:n+1) - f(1:n));
    
    u = un;

  endwhile

  err(n/100,1) = dx;
  err(n/100,2) = dx*sum(abs(u-uex));
  err(n/100,3) = sqrt(dx*sum((u-uex).^2));
  err(n/100,4) = max(abs(u-uex));

end

err

save "err_lw.dat" err




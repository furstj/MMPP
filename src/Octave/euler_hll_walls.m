clear all;

% Reseni smisene ulohy pro Eulerovy rovnice
%   rho_t + (rho u)_x = 0
%   (rho u)_t + (rho u^2 + p)_x = 0
%   (rho E)_t + [u (rho E + p)]_x = 0
%
%   p = (kappa-1) * rho * (E - u^2/2)
%
%   rovnice resime ve tvaru W_t + F(W)_x = 0, kde
%   W = [rho, m, e], F(W) = [rho u, rho u^2 + p, u (rho E + p)] 
%   kde m = rho u, e = rho E 
%
%   Na krajich intervalu predpokladam pevne steny
%
kappa = 1.4;

% Delka intervalu a pocet bunek
L = 1;
n = 200;

dx = L / n;
x = dx/2:dx:(L-dx/2);

W(1:3,1:n)  = 0;

% Pocatecni podminka
for i = 1:n
  if (x(i)<0.5)
    W(1:3,i) = [ 12; 0; 1e6/(kappa-1) ];
  else
    W(1:3,i) = [ 1.2; 0; 1e5/(kappa-1) ];
  end
end

plot(x, W(1,:),  "-b;rho;"); axis([0 1 0 20]);
disp("Stiskni enter pro pokracovani"); pause;
t = 0;

for iter = 1:(20*n)
    
  rho = W(1,:);
  u   = W(2,:) ./ rho;
  E   = W(3,:) ./ rho;

  p = (kappa-1) * rho .* (E - 0.5*u.*u);
  a = sqrt(kappa * p ./ rho);

  FF(1,:) = rho .* u;
  FF(2,:) = rho .* u.^2 + p;
  FF(3,:) = u .* (rho .* E + p);

  % Vypocet spektralniho polomeru jakobianu
  sigma = abs(u) + a;
  %
  
  dt = 0.4 * dx / max(sigma);
  
  % Numericky tok na leve a prave casti intervalu
  pb = p(1);
  F(:,1) = [0; pb; 0];

  pb = p(n);
  F(:,n+1) = [0; pb; 0];
  
  % Numericky tok f
  for i = 2:n
    sl = min( u(i-1)-a(i-1), u(i)-a(i) );
    sr = max( u(i-1)+a(i-1), u(i)+a(i) );

    Wl = W(:,i-1);
    Wr = W(:,i);
    
    Fl = FF(:,i-1);
    Fr = FF(:,i);
    
    if (0<=sl)
      F(:,i) = Fl;
    elseif (0<sr)
      F(:,i) = (sr*Fl - sl*Fr + sl*sr*(Wr-Wl)) / (sr-sl);
    else
      F(:,i) = Fr;
    end

  end
  
  W(:,1:n) = W(:,1:n) - dt/dx * (F(:,2:n+1) - F(:,1:n));
  t = t + dt;
  
  if (mod(iter,10)==0)
    a = sqrt(kappa*p./W(1,:));
    plot(x, W(1,:),  "-b;rho;", x, W(2,:)./W(1,:)./a, "-k;u/a;", x, p/1e5, "-r;p;"); axis([0 1 -2 15]);
    pause(0.2);
  end
end


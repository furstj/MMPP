clear all;

% Reseni pocatecni ulohy pro rovnice melke vody
%   h_t + (uh)_x = 0
%   (hu)_t + (hu^2 + gh^2/2)_x = 0
%
%   rovnice resime ve tvaru W_t + F(W)_x = 0, kde
%   W = [h, q], F(W) = [q, q^2/h + gh^2/2] kde q = hu
%

g = 10;

% Delka intervalu a pocet bunek
L = 1;
n = 400;

dx = L / n;
x = dx/2:dx:(L-dx/2);

W(1:2,1:n)  = 0;

% Pocatecni podminka
for i = 1:n
  if (x(i)<0.5)
    W(1,i) = 1.0;
  else
    W(1,i) = 0.01;
  end
end

plot(x, W(1,:)); axis([0 1 0 1.5]);
disp("Stiskni enter pro pokracovani"); pause;
t = 0;

for iter = 1:n
    
    
  % Vypocet "derivace" W
  DW(1:2,1) = 0;
  DW(1:2,n) = 0;
  for i = 2:n-1
    for k = 1:2
      if ( (W(k,i)-W(k,i-1))*(W(k,i+1)-W(k,i)) <= 0 )
	DW(k,i) = 0;
      else
	if ( abs(W(k,i)-W(k,i-1)) < abs(W(k,i+1)-W(k,i)) ) 
	  DW(k,i) = W(k,i)-W(k,i-1);
	else
	  DW(k,i) = W(k,i+1)-W(k,i);
	end
      end
    end
  end

  h = W(1,:);
  q = W(2,:);
  u = q ./ h;
  
  % Vypocet spektralniho polomeru jakobianu
  sigma = abs(u) + sqrt(g*h);
  %
  
  dt = 0.4 * dx / max(sigma);
  
  % Numericky tok na leve a prave casti intervalu
  F(:,1) = [q(1); q(1)^2/h(1) + g*h(1)^2/2] ;
  F(:,n+1) = [q(n); q(n)^2/h(n) + g*h(n)^2/2] ;
  
  % Numericky tok f
  for i = 2:n
    s = max(sigma(i-1),sigma(i));
    Wl = W(:,i-1) + DW(:,i-1)/2;
    Wr = W(:,i) - DW(:,i)/2;
    
    Fl = [Wl(2); Wl(2)^2/Wl(1) + g*Wl(1)^2/2];
    Fr = [Wr(2); Wr(2)^2/Wr(1) + g*Wr(1)^2/2];
    
    F(:,i) = (Fl+Fr)/2 - s/2*(Wr-Wl);
  end
  
  W(:,1:n) = W(:,1:n) - dt/dx * (F(:,2:n+1) - F(:,1:n));
  t = t + dt;
  
  if (mod(iter,10)==0)
    plot(x, W(1,:)); axis([0 1 0 1.5]);
    pause(1);
  end
end


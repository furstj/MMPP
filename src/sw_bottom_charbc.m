clear all;

% Reseni proudeni melke vody pres prekazku
%
%   h_t + (uh)_x = 0
%   (hu)_t + (hu^2 + gh^2/2)_x = -g h b'
%
%   rovnice resime ve tvaru W_t + F(W)_x = Q, kde
%   W = [h, q], F(W) = [q, q^2/h + gh^2/2] kde q = hu a 
%   Q = [0, -g h b' - cf*u*|u|]
%
%
%  Jakobian dF/dW = [ 0, 1; gh^2-u, 2u]
%  vlastni cisla  u +- sqrt(g*h)
%  vlastni vektory [1; u+-sqrt(g*h)]
%
%  charakteristicke promenne V = [u+2*c,u-2*c]     c =sqrt(g*h)

g = 10;
cf = 2;

% Rychlost vtoku m3/s
u0 = 0.5;
hL = 0.5;


% Delka intervalu a pocet bunek
L = 1;
n = 200;

dx = L / n;
x = dx/2:dx:(L-dx/2);

% tvar dna
for i = 1:n
    if (0.3 < x(i) && x(i)<0.7)  
      b(i) = (cos(pi*(x(i)-0.5)/0.2) + 1)/2 + (1-x(i))*0.2; 
      bx(i) = ( (cos(pi*(x(i)+dx/2-0.5)/0.2) + 1)/2 - 
		(cos(pi*(x(i)-dx/2-0.5)/0.2) + 1)/2 ) / dx - 0.2; 
    else
      b(i) = (1-x(i))*0.2;
      bx(i)= -0.2;
    end
end

W(1,1:n)  = max(hL,max(b)+0.1) - b;
W(2,1:n)  = W(1,:)*u0;

plot(x, b, x, W(1,:)+b); axis([0 1 0 3]);
disp("Stiskni enter pro pokracovani"); pause;
t = 0;

tend = 1;

Q(1:2,1:n) = 0;

for iter = 1:(40*n)
    
  h = W(1,:);
  q = W(2,:);
  u = q ./ h;
  
  % Vypocet spektralniho polomeru jakobianu
  sigma = abs(u) + sqrt(g*h);
  %
  
  dt = 0.4 * dx / max(sigma);

  
  % Numericky tok na vstupu (zadana rychlost)
  % u0 - 2c0 := u1 - 2c1
  c0 = (u0-u(1)) / 2.0 + sqrt(g*h(1));
  h0 = c0^2 / g;
  F(:,1) = [h0*u0; h0*u0^2 + g*h0^2/2] ;

  % Numericky tok na vystupu (konstantni hodnota na charakteristice)
  % uL + 2*cL := un + 2*cn 
  uL = u(n) + 2*(sqrt(g*h(n)) - sqrt(g*hL));
  F(:,n+1) = [hL*uL; hL*uL^2 + g*hL^2/2] ;
  
  % Numericky tok f
  for i = 2:n
    sl = min( u(i-1)-sqrt(g*h(i-1)), u(i)-sqrt(g*h(i)) );
    sr = min( u(i-1)+sqrt(g*h(i-1)), u(i)+sqrt(g*h(i)) );

    Wl = W(:,i-1);
    Wr = W(:,i);
    
    Fl = [Wl(2); Wl(2)^2/Wl(1) + g*Wl(1)^2/2];
    Fr = [Wr(2); Wr(2)^2/Wr(1) + g*Wr(1)^2/2];
    
    if (0<=sl)
      F(:,i) = Fl;
    elseif (0<sr)
      F(:,i) = (sr*Fl - sl*Fr + sl*sr*(Wr-Wl)) / (sr-sl);
    else
      F(:,i) = Fr;
    end

  end

  % Zdrojovy clen
  Q(2,1:n) = - g * h .* bx - cf*abs(u).*u;

  W(:,1:n) = W(:,1:n) - dt/dx * (F(:,2:n+1) - F(:,1:n)) + dt*Q;
  t = t + dt;
  
  if (mod(iter,10)==0)
    plot(x, b, x, W(1,:)+b); axis([0 1 0 3]);
    pause(0.1);
  end
end


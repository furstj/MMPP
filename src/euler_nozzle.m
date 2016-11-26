clear all;

% Reseni proudeni 1D dyzou
%   (A rho)_t + (A rho u)_x = 0
%   (A rho u)_t + ( A rho u^2 + A p)_x = p A_x
%   (A rho E)_t + [A u (rho E + p)]_x = 0
%
%   p = (kappa-1) * rho * (E - u^2/2)
%
%   rovnice resime ve tvaru A W_t + (A F(W))_x = Q, kde
%   W = [rho, m, e], F(W) = [rho u, rho u^2 + p, u (rho E + p)] 
%   kde m = rho u, e = rho E 
%
%   na vstupu zadvame celkovy tlak a celkovou teplotu
%   na vystupu zadavame staticky tlak
%
kappa = 1.4;

% Delka intervalu a pocet bunek
L = 1;
n = 200;

dx = L / n;
x = dx/2:dx:(L-dx/2);

% Prurezy dyzy na hranicich bunek  A(x) = 1 + (x-0.5)^2
A(1:n+1) = ((0:dx:L)-L/2) .^2 + 1;
Ai(1:n) = (A(1:n) + A(2:n+1)) / 2;

W(1:3,1:n) = 0;
Q(1:3,1:n) = 0;

% Pocatecni podminka
for i = 1:n
    W(1:3,i) = [ 1.1; 0; 1e5/(kappa-1) ];
end

plot(x, W(1,:)); axis([0 1 0 1.2]);
disp("Stiskni enter pro pokracovani"); pause;
t = 0;

Rezidua = [];

for iter = 1:20000
    
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
  pTot = 1.e5;
  rhoTot = 1.2;

  pb =  min(p(1),pTot);
  M2 = ( (pb/pTot) ^ ((1-kappa)/kappa) - 1 ) * 2/(kappa-1);
  rhob = rhoTot * (1+(kappa-1)/2*M2) ^ (1/(1-kappa));
  cb = sqrt(kappa*pb/rhob);
  ub = cb*sqrt(M2);
  rEb = pb/(kappa-1) + 0.5*rhob*ub^2;
  F(:,1) = [rhob*ub; rhob*ub^2+pb; ub*(rEb + pb)];

  p2 = 0.75 * pTot;

  pb = p2;
  rhob = rho(n);
  ub = u(n);
  rEb = pb/(kappa-1) + 0.5*rhob*ub^2;
  F(:,n+1) = [rhob*ub; rhob*ub^2+pb; ub*(rEb + pb)];
  
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
  Q(2,:) = p .* (A(2:n+1)-A(1:n))/dx;


  for k=1:3
    F(k,:) = A .* F(k,:);
    Wn(k,1:n) = W(k,1:n) - dt/dx./Ai .* (F(k,2:n+1) - F(k,1:n)) + dt./Ai.*Q(k,:); 
  end
  
  if (mod(iter,100)==0)
    Rezidua = [ Rezidua,  sqrt(sum((Wn-W).**2,dim=2)) / dt];
    plot(x, p/1e5); axis([0 1 0 1.2]);
    pause(0.01);
  end

  W = Wn;
  t = t + dt;
  
end


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

for iter = 1:(5*n)

    h = W(1,:);
    q = W(2,:);
    u = q ./ h;

    % Vypocet fyzikalniho toku
    FF(1,:) = q;
    FF(2,:) = q.^2./h + g*h.^2/2;

    % Vypocet spektralniho polomeru jakobianu
    sigma = abs(u) + sqrt(g*h);
    %

    dt = 0.4 * dx / max(sigma);

    % Numericky tok na leve a prave casti intervalu
    W0  = [ h(1); -q(1) ];
    FF0 = [-FF(1,1); FF(2,1)] ;
    s = sigma(1);
    F(:,1) = (FF0+FF(:,1))/2 - s/2*(W(:,1)-W0);

    WL  = [ h(n); -q(n) ];
    FFL = [-FF(1,n); FF(2,n)] ;
    s = sigma(n);
    F(:,n+1) = (FFL+FF(:,n))/2 - s/2*(WL-W(:,n));

    % Numericky tok f
    for i = 2:n
	s = max(sigma(i-1),sigma(i));
	F(:,i) = (FF(:,i-1)+FF(:,i))/2 - s/2*(W(:,i)-W(:,i-1));
    end

    W(:,1:n) = W(:,1:n) - dt/dx * (F(:,2:n+1) - F(:,1:n));
    t = t + dt;

    if (mod(iter,10)==0)
      plot(x, W(1,:)); axis([0 1 0 1.5]);
      pause(0.1);
    end
end




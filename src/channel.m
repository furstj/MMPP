% Simulace proudeni nevazke tekutiny v kanale s promennym prurezem

clear all;

% Velikost site (pocty bunek)
ni = 300;
nj = 100;

% Okrajove podminky (+nulovy uhel na vstupu) a parametry plynu
pTot = 1.e5;
rhoTot = 1.2;
p2 = 0.737*pTot;
kappa = 1.4;
 
% Sit
x(:,1) = 0:3/ni:3;
x(:,nj+1) = x(:,1);

for i=1:ni+1
    if (x(i,1) < 1)
       y(i,1) = 0;
    elseif (x(i,1) < 2)
      h=0.1;
      r = (0.25+h^2)/(2*h);
      y(i,1) = sqrt( r^2-(x(i)-1.5)^2) + h - r;
    else
      y(i,1) = 0;
    end
end
y(:,nj+1)=1;

for j=2:nj
    x(:,j) = x(:,1);
    y(:,j) = ( (j-1)*y(:,nj+1) + (nj+1-j)*y(:,1) ) / nj;
end

for i=1:ni
    for j=1:nj
      Vol(i,j) = 0.5*cross(
			 [x(i+1,j+1)-x(i,j);y(i+1,j+1)-y(i,j);0], 
			 [x(i,j+1)-x(i+1,j);y(i,j+1)-y(i+1,j);0])(3);
    end
end
xc = (x(1:ni,1:nj)+x(2:ni+1,1:nj)+x(2:ni+1,2:nj+1)+x(1:ni,2:nj+1)) / 4.0;
yc = (y(1:ni,1:nj)+y(2:ni+1,1:nj)+y(2:ni+1,2:nj+1)+y(1:ni,2:nj+1)) / 4.0;

pcolor(x,y, ones(ni+1,nj+1)); axis([0 3 -0.7 1.4 0 1]); 
disp("Obrazek site, stiskni enter pro pokracovani")
%pause

% Pocatecni podminka
W(1:ni,1:nj,1) = rhoTot;
W(1:ni,1:nj,2) = 0;
W(1:ni,1:nj,3) = 0;
W(1:ni,1:nj,4) = pTot/(kappa-1);

t = 0;
Rez = [];
for iter = 1:100000
    
  rho  = W(:,:,1);
  u    = W(:,:,2) ./ rho;
  v    = W(:,:,3) ./ rho;
  rhoE = W(:,:,4);
  
  p = (kappa-1) * (rhoE - 0.5*rho.*(u.^2+v.^2));
  
  a = sqrt(kappa*p./rho);
  
  Sx = y(2:ni,2:nj+1) - y(2:ni,1:nj); 
  Sy = x(2:ni,1:nj) - x(2:ni,2:nj+1); 
  magS = sqrt(Sx.*Sx + Sy.*Sy);
  
  ql = u(1:ni-1,:) .* Sx + v(1:ni-1,:) .* Sy;
  qr = u(2:ni,:) .* Sx + v(2:ni,:) .* Sy;
  
  sl = min(ql - magS.*a(1:ni-1,:), qr - magS.*a(2:ni,:));
  sr = min(ql + magS.*a(1:ni-1,:), qr + magS.*a(2:ni,:));
  
  sl = min(sl,0);
  sr = max(sr,0);
  
  F1l(:,:,1) = W(1:ni-1,:,1) .* ql;
  F1l(:,:,2) = W(1:ni-1,:,2) .* ql + p(1:ni-1,:).*Sx;
  F1l(:,:,3) = W(1:ni-1,:,3) .* ql + p(1:ni-1,:).*Sy;
  F1l(:,:,4) = (W(1:ni-1,:,4) + p(1:ni-1,:)) .* ql;
  
  F1r(:,:,1) = W(2:ni,:,1) .* qr;
  F1r(:,:,2) = W(2:ni,:,2) .* qr + p(2:ni,:).*Sx;
  F1r(:,:,3) = W(2:ni,:,3) .* qr + p(2:ni,:).*Sy;
  F1r(:,:,4) = (W(2:ni,:,4) + p(2:ni,:)) .* qr;
  
  for k=1:4
    F1(2:ni,:,k) = (sr.*F1l(:,:,k) - sl.*F1r(:,:,k) 
		    + sl.*sr.*(W(2:ni,:,k)-W(1:ni-1,:,k))) ./ (sr-sl);
  end

  % Vypocet na vstupni hranici (musi byt svisla!)
  for j = 1:nj
    i = 1;
    S = [ y(i,j+1)-y(i,j), x(i,j)-x(i,j+1) ];
    magS = sqrt(sum(S.^2));
      
    pb =  min(p(1,j),pTot);
    M2 = ( (pb/pTot) ^ ((1-kappa)/kappa) - 1 ) * 2/(kappa-1);
    rhob = rhoTot * (1+(kappa-1)/2*M2) ^ (1/(1-kappa));
    cb = sqrt(kappa*pb/rhob);
    ub = cb*sqrt(M2);
    rEb = pb/(kappa-1) + 0.5*rhob*ub^2;
    F1(1,j,:) = [rhob*ub; rhob*ub^2+pb; 0; ub*(rEb + pb)] * magS;
    
    % Vypocet na vystupni hranici (musi byt svisla!)
    i = ni+1;
    S = [ y(i,j+1)-y(i,j), x(i,j)-x(i,j+1) ];
    magS = sqrt(sum(S.^2));
    pb = p2;
    rhob = rho(ni,j);
    ub = u(ni,j);
    vb = v(ni,j);
    rEb = pb/(kappa-1) + 0.5*rhob*(ub^2+vb^2);
    F1(ni+1,j,:) = [rhob*ub; rhob*ub^2+pb; rhob*ub*vb; ub*(rEb + pb)]*magS;
  end

  
  % Vypocet toku pres hranice j-1/2  F2(i,j) = flux(W(:,i,j+1),W(:,i,j)) 
  Sx = y(1:ni,2:nj) - y(2:ni+1,2:nj); 
  Sy = x(2:ni+1,2:nj) - x(1:ni,2:nj); 
  magS = sqrt(Sx.*Sx + Sy.*Sy);
  
  ql = u(:,1:nj-1) .* Sx + v(:,1:nj-1) .* Sy;
  qr = u(:,2:nj) .* Sx + v(:,2:nj) .* Sy;
  
  sl = min(ql - magS.*a(:,1:nj-1), qr - magS.*a(:,2:nj));
  sr = min(ql + magS.*a(:,1:nj-1), qr + magS.*a(:,2:nj));
  
  sl = min(sl,0);
  sr = max(sr,0);
  
  F2l(:,:,1) = W(:,1:nj-1,1) .* ql;
  F2l(:,:,2) = W(:,1:nj-1,2) .* ql + p(:,1:nj-1).*Sx;
  F2l(:,:,3) = W(:,1:nj-1,3) .* ql + p(:,1:nj-1).*Sy;
  F2l(:,:,4) = (W(:,1:nj-1,4) + p(:,1:nj-1)) .* ql;
  
  F2r(:,:,1) = W(:,2:nj,1) .* qr;
  F2r(:,:,2) = W(:,2:nj,2) .* qr + p(:,2:nj).*Sx;
  F2r(:,:,3) = W(:,2:nj,3) .* qr + p(:,2:nj).*Sy;
  F2r(:,:,4) = (W(:,2:nj,4) + p(:,2:nj)) .* qr;
  
  for k=1:4
    F2(:,2:nj,k) = (sr.*F2l(:,:,k) - sl.*F2r(:,:,k) 
		    + sl.*sr.*(W(:,2:nj,k)-W(:,1:nj-1,k))) ./ (sr-sl);
  end


  for i = 1:ni
      % pevne steny
      j = 1;
      S = [ y(i,j)-y(i+1,j), x(i+1,j)-x(i,1) ];
      magS = sqrt(sum(S.^2));
      F2(i,j,:) = [0; S(1); S(2); 0] * p(i,1);

      j = nj+1;
      S = [ y(i,j)-y(i+1,j), x(i+1,j)-x(i,1) ];
      magS = sqrt(sum(S.^2));
      F2(i,j,:) = [0; S(1); S(2); 0] * p(i,j-1);

    end


    % Zjednoduseny odhad casoveho kroku
    dx = 3 / ni;
    dy = 3 / nj;
    dt = 0.4 / ( max(max(abs(u)+a))/dx + max(max(abs(v)+a))/dy);

    for k = 1:4
	Wn(:,:,k) = W(:,:,k) - dt./Vol .* (F1(2:ni+1,:,k)-F1(1:ni,:,k) + 
					   F2(:,2:nj+1,k)-F2(:,1:nj,k));
    end

    if ( mod(iter,10)==0)
      pcolor(xc,yc,p); axis([0 3 -0.7 1.4 0 1]); shading interp;       
      pause(0.001);
    end

    thisRez = sqrt(sum(sum((Wn-W).^2)))/dt/ni/nj;

    Rez = [Rez; [thisRez(1,1,1),thisRez(1,1,2),thisRez(1,1,3),thisRez(1,1,4)]];
    t = t + dt;
    W = Wn;

end


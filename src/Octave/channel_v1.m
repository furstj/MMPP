% Simulace proudeni nevazke tekutiny v kanale s promennym prurezem

clear all;

% Velikost site (pocty bunek)
ni = 60;
nj = 20;

% Okrajove podminky (+nulovy uhel na vstupu) a parametry plynu
pTot = 1.e5;
rhoTot = 1.2;
p2 = 0.675*pTot;
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
      Vol(1,i,j) = 0.5*cross(
			 [x(i+1,j+1)-x(i,j);y(i+1,j+1)-y(i,j);0], 
			 [x(i,j+1)-x(i+1,j);y(i,j+1)-y(i+1,j);0])(3);
    end
end
xc = (x(1:ni,1:nj)+x(2:ni+1,1:nj)+x(2:ni+1,2:nj+1)+x(1:ni,2:nj+1)) / 4.0;
yc = (y(1:ni,1:nj)+y(2:ni+1,1:nj)+y(2:ni+1,2:nj+1)+y(1:ni,2:nj+1)) / 4.0;

pcolor(x,y, ones(ni+1,nj+1)); axis([0 3 -0.7 1.4 0 1]); 
disp("Obrazek site, stiskni enter pro pokracovani")
pause

% Pocatecni podminka
W(1,1:ni,1:nj) = rhoTot;
W(2,1:ni,1:nj) = 0;
W(3,1:ni,1:nj) = 0;
W(4,1:ni,1:nj) = pTot/(kappa-1);

t = 0;
for iter = 1:1000

    rho  = reshape(W(1,:,:), [ni,nj]);
    u    = reshape(W(2,:,:),[ni,nj]) ./ rho;
    v    = reshape(W(3,:,:),[ni,nj]) ./ rho;
    rhoE = reshape(W(4,:,:), [ni,nj]);

    p = (kappa-1) * (rhoE - 0.5*rho.*(u.^2+v.^2));

    a = sqrt(kappa*p./rho);

    % Vypocet toku pres hranice i-1/2  F1(i,j) = flux(W(:,i-1,j),W(:,i,j)) 
    for j = 1:nj
	for i = 2:ni
	    S = [ y(i,j+1)-y(i,j), x(i,j)-x(i,j+1) ];
	    magS = sqrt(sum(S.^2));

	    Wl = W(:,i-1,j);
	    ql = u(i-1,j)*S(1) + v(i-1,j)*S(2);
	    F1l = ql * Wl + p(i-1,j)*[0; S(1); S(2); ql];

	    Wr = W(:,i,j);
	    qr = u(i,j)*S(1) + v(i,j)*S(2);
	    F1r = qr * Wr + p(i,j)*[0; S(1); S(2); qr];

	    sl = min(ql - magS*a(i-1,j), qr - magS*a(i,j));
	    sr = max(ql + magS*a(i-1,j), qr + magS*a(i,j));

	    if (0<=sl)
	      F1(:,i,j) = F1l;
	    elseif (0<sr)
	      F1(:,i,j) = (sr*F1l - sl*F1r + sl*sr*(Wr-Wl)) / (sr-sl);
	    else
	      F1(:,i,j) = F1r;
	    end
	end

	% Vypocet na vstupni hranici (musi byt svisla!)
	i = 1;
	S = [ y(i,j+1)-y(i,j), x(i,j)-x(i,j+1) ];
	magS = sqrt(sum(S.^2));

	pb =  min(p(1,j),pTot);
	M2 = ( (pb/pTot) ^ ((1-kappa)/kappa) - 1 ) * 2/(kappa-1);
	rhob = rhoTot * (1+(kappa-1)/2*M2) ^ (1/(1-kappa));
	cb = sqrt(kappa*pb/rhob);
	ub = cb*sqrt(M2);
	rEb = pb/(kappa-1) + 0.5*rhob*ub^2;
	F1(:,1,j) = [rhob*ub; rhob*ub^2+pb; 0; ub*(rEb + pb)] * magS;

	% Vypocet na vystupni hranici (musi byt svisla!)
	i = ni+1;
	S = [ y(i,j+1)-y(i,j), x(i,j)-x(i,j+1) ];
	magS = sqrt(sum(S.^2));
	pb = p2;
	rhob = rho(ni,j);
	ub = u(ni,j);
	vb = v(ni,j);
	rEb = pb/(kappa-1) + 0.5*rhob*(ub^2+vb^2);
	F1(:,ni+1,j) = [rhob*ub; rhob*ub^2+pb; rhob*ub*vb; ub*(rEb + pb)]*magS;
    end

    % Vypocet toku pres hranice j-1/2  F2(i,j) = flux(W(:,i,j+1),W(:,i,j)) 
    for i = 1:ni
      for j = 2:nj
	S = [ y(i,j)-y(i+1,j), x(i+1,j)-x(i,1) ];
	magS = sqrt(sum(S.^2));
	
	Wl = W(:,i,j-1);
	ql = u(i,j-1)*S(1) + v(i,j-1)*S(2);
	F2l = ql * Wl + p(i,j-1)*[0; S(1); S(2); ql];
	
	Wr = W(:,i,j);
	qr = u(i,j)*S(1) + v(i,j)*S(2);
	F2r = qr * Wr + p(i,j)*[0; S(1); S(2); qr];
	
	sl = min(ql - magS*a(i,j-1), qr - magS*a(i,j));
	sr = max(ql + magS*a(i,j-1), qr + magS*a(i,j));
	
	if (0<=sl)
	  F2(:,i,j) = F2l;
	elseif (0<sr)
	  F2(:,i,j) = (sr*F2l - sl*F2r + sl*sr*(Wr-Wl)) / (sr-sl);
	else
	  F2(:,i,j) = F2r;
	end
      end
      
      % pevne steny
      j = 1;
      S = [ y(i,j)-y(i+1,j), x(i+1,j)-x(i,1) ];
      magS = sqrt(sum(S.^2));
      F2(:,i,j) = [0; S(1); S(2); 0] * p(i,1);

      j = nj+1;
      S = [ y(i,j)-y(i+1,j), x(i+1,j)-x(i,1) ];
      magS = sqrt(sum(S.^2));
      F2(:,i,j) = [0; S(1); S(2); 0] * p(i,j-1);

    end


    % Zjednoduseny odhad casoveho kroku
    dx = 3 / ni;
    dy = 3 / nj;
    dt = 0.4 / ( max(max(abs(u)+a))/dx + max(max(abs(v)+a))/dy);

    for k = 1:4
	Wn(k,:,:) = W(k,:,:) - dt./Vol .* (F1(k,2:ni+1,:)-F1(k,1:ni,:) + 
					   F2(k,:,2:nj+1)-F2(k,:,1:nj));
    end

    if ( mod(iter,10)==0)
      pcolor(xc,yc,p); axis([0 3 -0.7 1.4 0 1]);        
      pause(0.01);
    end

    t = t + dt;
    W = Wn;

end


% Quick test to see how Murrays model works for other mass ratios
clear

syms m1 k b w Q E t

m1 = 0.5;

m2 = 1 - m1;

k1 = k*(m1+m2)/m2;

b1 = b*(m1+m2)/m2;

y1 = -b1/(-2*m1); % Eqn 12, soln to quaratic

A1 = Q*E/(b1*1i*w); %Eqn 16, at resonance

v1 = A1*1i*w;

P1 = (Q*E)^2/(2*y1*m1) % At resonance

%%% Note full solution in other


k2 = k*(m1+m2)/m1;

b2 = b*(m1+m2)/m1;

y2 = -b2/(-2*m2); % Eqn 12, soln to quaratic

A2 = Q*E/(b2*1i*w); %Eqn 16, at resonance

v2 = A2*1i*w;

P2 = (Q*E)^2/(2*y2*m2)

P = P1 + P2



mr = (m1*m2)/(m1+m2);

y = -b/(-2*mr); % Eqn 12, soln to quaratic

A = Q*E/(b*1i*w); %Eqn 16, at resonance

v = A*1i*w;

P = (Q*E)^2/(2*y*mr)


% From Sun

function P = engine(t, C_star, rho, Dt, w) 
% response surface of P(t) in the 
% solid rocket engine problem

func = {
%     @(x1,x2,x3) x1.^2;
%     @(x1,x2,x3) x2.^2;
    @(x1,x2,x3) x3.^2;
%     @(x1,x2,x3) x1.*x2;
    @(x1,x2,x3) x1.*x3;
    @(x1,x2,x3) x2.*x3;
%     @(x1,x2,x3) x1;
%     @(x1,x2,x3) x2;
    @(x1,x2,x3) x3;
    @(x1,x2,x3) 1; };

w_t = 0;
for jj = 1:length(func)
    w_t = w_t + w(jj,:) * func{jj}(C_star, rho, Dt);
end
P = w_t(1).*t.^2 + w_t(2).*t + w_t(3);
end
function t=findthresh(d,l,b,deltax,p)

myfun = @(t,d,l,b,deltax,p) solve_pde_1d_prob(t,d,l,b,deltax)-p;  % parameterized function

fun = @(x) myfun(x,d,l,b,deltax,p);
if sign(fun(.1))~=sign(fun(10))
    t = fzero(fun,[.1 10]);
else 
    t=fzero(fun,2);
    v=[.1 t 10];
    q=[fun(.1) fun(t) fun(10)];
    [~,m]=min(abs(q));
    t=v(m);
end

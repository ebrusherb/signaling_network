function new = transform(p)
    x=p(:,1);
    %y=p(:,2);
    %z=p(:,3);
    y=p(:,3);
    z=p(:,2);
    xprime=1-x+z;
    yprime=sqrt(3)*y;
    new=[xprime yprime];
end
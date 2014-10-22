Ni=4
Nj=5
u=array(,(Ni-1)*(Nj-1))
A=array(,c(Ni-1)*(Nj-1),(Ni-1)*(Nj-1))

fval=0
f=array(fval,(Ni-1)*(Nj-1))

deltax=.01
Lx=0
Ux=10
xvals=seq(Lx,Ux,by=deltax)

deltay=.01
L

bl=1
br=0
bd=0
du=1
b=array(0,(Ni-1)*(Nj-1))
down=1+(1:(Ni-2))*(Ni-1)
up=Ni-1+(1:(Ni-2))*(Ni-1)
left=2:(Ni-2)
right=(Nj-2)*(Ni-1)+2:(Ni-2)


library(colorspace)
setwd("/Users/eleanorbrush/Dropbox/signaling_network")

n=101;

cols=sequential_hcl(n=n,h=120,l=c(30,100),c.=c(100,0),power=.5)
cols=col2rgb(cols)
cols=cols/255
cols=matrix(cols,nrow=3)
cols=cols[,n:1]
write(cols,file='green.txt',ncolumns=3)
 
cols=sequential_hcl(n=n,h=10,l=c(30,100),c.=c(100,0),power=.5)
cols=col2rgb(cols)
cols=cols/255
cols=matrix(cols,nrow=3)
cols=cols[,n:1]
write(cols,file='red.txt',ncolumns=3)
 
cols=sequential_hcl(n=n,h=260,l=c(30,100),c.=c(100,0),power=.5)
cols=col2rgb(cols)
cols=cols/255
cols=matrix(cols,nrow=3)
cols=cols[,n:1]
write(cols,file='blue.txt',ncolumns=3)
 
cols=sequential_hcl(n=n,l=c(30,100),c.=c(0,0),power=.5)
cols=col2rgb(cols)
cols=cols/255
cols=matrix(cols,nrow=3)
cols=cols[,n:1]
write(cols,file='black.txt',ncolumns=3)


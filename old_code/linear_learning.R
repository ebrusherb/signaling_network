binsig<-function(name){ #puts 1 where there nonzero entries of name and 0s elsewhere
	n<-dim(name)[1]
	x<-array(0,c(n,n))
for(i in 1:n){for(j in 1:n){if(name[i,j]!=0){x[i,j]=1}}}
return(x)
}

laplace<-function(i,mat){
	M=mat
	diag(M)=1
	w=c(which(M[i,]!=0))
	M=M[w,w]
	B=binsig(M)
	new=array(0,c(length(w)-1,length(w)-1))
	for(j in 1:length(w)-1){
		new[j,j]=-sum(B[j,])-B[j,length(w)]
		new[j,-j]=B[j,-c(j,length(w))]-B[j,length(w)]
		}
	diag(new)=diag(new)+1
	return(new)
	}

mat=matrix(sample(0:10,100,replace=TRUE),nrow=8)	
N=50
mat=matrix(0,nrow=N,nc=N)
mat[sample(1:N^2,floor(N^2/2))]=sample(1:N^2,floor(N^2/2),replace=TRUE)
mat[lower.tri(mat)]=t(mat)[lower.tri(mat)]

mat=prox

mat=get.adjacency(erdos.renyi.game(100,p=.4,type="gnp"))

N=dim(mat)[2]
learning<-array(,N)
for(i in 1:N){
	M=laplace(i,mat)
	values<-Re(eigen(M)$values)
	learning[i]=max(values)
	}
	
L=laplace(1,prox)
eigen(L)$values
x=rep(0,31)
y=runif(31,0,0.001)
for(i in 1:10){
	#print(sum((L%*%y)^2))
	print(y[1])
	y=L%*%y
	}
	
lap<-function(i,mat){
	B=binsig(mat)
	diag(B)=-RR(B)
	return(B)
	}
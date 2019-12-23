#set.seed(1)
stoch_cusp=function(x,b,a,sdW,dt)
{
#  dt=.1
  x=x-dt*(x^3-b*x-a)+rnorm(1,0,sdW)
  x[is.nan(x)]=0
  x
}

stoch_cusp_run=function(x0,b,a,s,n,dt)
{
o=rep(0,n)
o[1]=x0
for(i in 2:n)
  o[i]=stoch_cusp(o[i-1],b,a,s,dt)
o
}

pdf("figure2.pdf",h=6,w=12)
layout(matrix(1:2,1,2))
set.seed(10)
dt=.15
s=.01
a=0;b=1
n=2000
x0=0
I_min=-.5
# illustion of stochastic cusp
plot(stoch_cusp_run(x0,b+I_min,a,s,n,dt)[-1:-50],type='l',bty='n',ylab='Opinion',xlab="time",,ylim=c(-1.5,1.5))

n=400;dt=.2;bb=c(1.5,1,.5,0)
a=seq(-.5,.5,le=n)
a=c(a,rev(a))
bi=0
so=.01
for(b in bb)
{
  bi=bi+1
  o=rep(0,n)
  runs=1000
  for(i in 2:(2*n))
    o[i]=stoch_cusp_run(o[i-1],b+I_min,a[i],s=so,n=runs,dt)[runs] # o[i-1]-dt*(o[i-1]^3-b*o[i-1]-a[i])
  if(b==bb[1])plot(o[-1:-10]~a[-1:-10],xlim=c(-.65,.65),type='l',bty='n',ylab='Opinion',xlab="Information (I)",main="") else lines(o[-1:-10]~a[-1:-10],col=bi,lty=bi)
}
legend("bottomright",lty=1:4,col=1:4,legend=paste("A = ",bb),cex = 0.9,bty='n')

dev.off()



set.seed(1)
pdf('figure3_involvement.pdf',h=4,w=8)
"Dynamics of involvement"
layout(1)
n=500

b=.9
x=rep(0,n)
x[1]=.01
p=1/(100^2)
p=.01
d=.01
I_equi=.5
a=I_equi*d/p
#b=(1-(a))/2

u=sample(0:1,n,T,prob=c(1-p,p))
for(i in 2:n)
  x[i]=(1-d)*x[i-1]+a*u[i]
plot(x,type='l',bty='n',main="",xlab='time',ylab="Attention")
lines(u/20,col='red')

dev.off()

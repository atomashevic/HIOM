require(igraph)
require(colorspace)
require(svMisc)
library(ggplot2)
library(Hmisc)
library(plotrix)
library(diptest)
library(png)
library(qgraph)
library(Rmisc) 
library(tidyverse)


hartigan_d= function(x)
{
  dip=dip.test(x)
  dip.p=dip$p.value
  dip=dip$statistic
  dipstar=4-as.numeric(cut(dip.p,c(-Inf,.001,.01,.05,Inf)))
  if(dipstar!=0) dipstar=paste(rep('*',dipstar),sep='',collapse="") else dipstar=""
  return(list(dip,dipstar))
}

std.error <- function(x) sd(x)/sqrt(length(x)) 

sample = function(x, size, replace = F, prob = NULL) {
  if (length(x) == 1) return(x)
  base::sample(x, size = size, replace = replace, prob = prob)
}

information_update=function(i1,i2,a1,a2,persuasion,r_min,info_update=T) # compute weighted man of information
{
  r=r_min+(1-r_min)/(1+exp(-1*persuasion*(a1-a2))) # resistance
  inf=r*i1+(1-r)*i2
  if(!info_update) inf=i1
  return(inf)
}

stoch_cusp=function(n,x,b,a,s_O,maxwell_convention=F) # updates opinion according to cusp
{
  dt=.01
  x=x-dt*(x^3-b*x-a)+rnorm(n,0,s_O)
  if(maxwell_convention) x[x*a<0]=-x[x*a<0]
  x[is.nan(x)]=0
  x
}

make_network = function(model='SBM',n=400,clusters=10,p_within=.01,p_between=.001,rewiring=.02)
{
  if(model=='WS')
  {
    g <- gg <- sample_smallworld(1, n,  1,rewiring);  # The Watts-Strogatz small-world model
    l <-layout_nicely(g)
    adj <- get.adjacency(g)
  }
  
  if(model=='SBM')
  {
    pm=matrix(p_between,clusters,clusters);
    diag(pm)=p_within
    g <- gg <- sample_sbm(n, pref.matrix=pm, block.sizes=rep(n/clusters,clusters))
    l <-layout_nicely(g)
    adj <- get.adjacency(g)
  }
  
  if(model=='lattice')
  {
    g=gg=make_lattice(c(sqrt(n), sqrt(n), 1),nei=1);
    l = layout_on_grid(g);
    adj <- get.adjacency(g)
  }
  
  if(model=='SBM_4')
  {  
    pm=matrix(p_between,clusters,clusters);
    diag(pm)=p_within
    g <- gg <- sample_sbm(.5*n, pref.matrix=pm, block.sizes=rep((.5*n)/clusters,clusters))
    l <-layout_nicely(g)
    g=g %du% g
    l1=l
    l1[,1]=l1[,1]
    l1[,2]=l1[,2]+45 
    l=rbind(l,l1)
    adj <- get.adjacency(g)
  }
  
  return(list(g,l,adj,gg))
}


plot.graph = function(adj,l,opin,inform,atten,CD,title="",shape='circle',c.title=1.5)
{
  vertex.c=diverging_hsv(256,s=2,v=1)[as.numeric(cut(opin,seq(min.o,max.o,length=255)))]
  vertex.border.c=rep('red',length(opin))
  vertex.border.c[inform<0]='blue'
  vertex.border.c[opin*inform>=0 | CD>0]='black'
  border.w=rep(1,length(opin))
  border.w[opin*inform>=0| CD>0]=.5
  par(mar=c(0,0,0,0))
  
  qgraph(adj, layout = l, vsize = 0.4 * (2+3*atten), labels = FALSE, color=vertex.c,
         shape = shape, esize = 0.5, edge.color = rgb(0.5,0.5,0.5, 0.5),border.color=vertex.border.c,
         border.width = border.w, plot = T, rescale = T,title=title,title.cex=c.title)
  par(mar=c(5,4,4,2));
}

plot.histo= function(x,min,max,xlab='')
{
  hist.c=rev(diverging_hsv(40,s=2,v=1))
  par(mar=c(3,0,0,0));par(mgp = c(1, 0, 0))
  hist(x,main='',xlab=xlab,ylab='',breaks=seq(min,max,le=40),col=hist.c,cex.lab=1.1,axes=F);
  axis(1,at=0,labels = '')
  par(mar=c(5,4,4,2));par(mgp = c(3, 1, 0))
}

##############

scenario=1

pdfplot=F # if TRUE to plot to pdf files

PNG=T # save png for gif
if(PNG) unlink(paste0("pngplots_",scenario,"/*"))

#layout for plots
if(scenario %in% c(1)) layout.m=layout.basis=layout(1)
if(scenario %in% c(2,3,31)) layout.m=layout.basis=matrix(c(rep(1,8),0,2,2,3,3,4,4,0),2,8,byrow=T)
if(scenario %in% c(4,41)) layout.m=layout.basis=matrix(c(rep(1,16),2,3,4,5),4,5)
if(scenario %in% c(1)) heights=heights.basis=1
if(scenario %in% c(2,3,31)) heights=heights.basis=c(4,1)
if(scenario %in% c(4,41)) heights=heights.basis=c(1,1)

if(scenario>30) simulation=T else simulation=F 

if(simulation) {plotting=F} else{plotting=T;sim_values=sim_values2=1;sim_values3='SBM'}  

if(scenario==31) {sim_values=rep(seq(0,.4,by=.05),each=5);sim_var=c('d_A');
sim_values2=rep(seq(.0,1,by=.1),each=5);sim_var2=c('r_min')
sim_values3=c(1,2);sim_var3=c('network')}
if(scenario==41) {sim_values=rep(seq(0,.002,by=.0005),each=5);sim_var=c('p(pertubation)')
sim_values2=rep(seq(0,.5,by=.1),each=5);sim_var2=c('t_d')
sim_values3=c(1,2);sim_var3=c('network')}

datasim=matrix(NA,length(sim_values)*length(sim_values2)*length(sim_values3),7)
sim_i=0 #counter in simulation loop

for(sim_value in sim_values)
{
  for(sim_value2 in sim_values2)
  {
    for(sim_value3 in sim_values3)
    {
      sim_i=sim_i+1
      print(sim_i)
      
      if(scenario==1) #the cusp of cusp; removed from manuscript
      {
        Ni= 5000 # max iterations
        
        if(pdfplot)
        { 
          set.seed(1);
          pdfname=paste('CAcusp_rev_s',scenario,'.pdf',collapse='',sep='')
          pdf(pdfname,h=8,w=11);
          plot_iteration=c(1,(1/10)*Ni,Ni)
          layout.m=matrix(1:3,3,1)
        }
        
        
        
        n1=200;n2=200  # dimensions of lattice
        n=n1*n2
        lattice=T
        g=make_lattice(c(n1, n2, 1),nei=1);l = layout_on_grid(g); # lattice
        
        sd_noise_information=.0
        info_update=F # if false information is not updated
        persuasion=1
        r_min=0
        
        s_O=.00001
        maxwell_convention=F
        
        attention_star=1
        min_attention=-.5 # to include continuous change in O as function of K
        delta_attention=0 # increase in attention when interacting
        decay_attention=delta_attention/(attention_star*(n^2))
        deffuant_c=Inf
        
        # initial values for A and I
        k_init=.2
        i_init_min=.5
        i_init_max=1
        information=rep(seq(-1*k_init,k_init,length=n1),each=n2)
        attention=rep(seq(i_init_min,i_init_max,length=n1),n2)
        opinion=rnorm(n,0,.3)
      }
      
      if(scenario==2)
      {
        Ni=15000# 5000 # max iterations
        if(pdfplot)
        { 
          set.seed(1);
          pdfname=paste('CAcusp_rev_s',scenario,'.pdf',collapse='',sep='')
          pdf(pdfname,h=8,w=11);
          layout.m=cbind(layout.basis,layout.basis+4*(layout.basis>0));layout.m=rbind(layout.m,layout.m+8*(layout.m>0))
          heights=c(3,1,3,1)
          plot_iteration=c(1,(1/3)*Ni,(2/3)*Ni,Ni)
        }
        
        n=400
        network=make_network(model='WS',n=n,clusters=10,p_within=.2,p_between=.001,rewiring=.02)
        g=network[[1]];l=network[[2]];adj=network[[3]]
        
        sd_noise_information=.0005
        info_update=T
        persuasion=1
        r_min=0
        
        s_O=.01
        maxwell_convention=F
        
        attention_star=1
        min_attention=-.5
        delta_attention=0
        decay_attention=delta_attention/(attention_star*(n^2))
        deffuant_c=Inf
        
        # all slighty positive low attention
        information = rnorm(n, 0, .3)
        attention = runif(n, 1, 1)
        opinion=rnorm(n,0,.2)
        
        # setup pdf files
      }
      
      if(scenario==3)
      {
        Ni= 25000# 20000 # max iterations
        
        if(pdfplot) 
        {
          set.seed(1);
          pdfname=paste('CAcusp_rev_s',scenario,'.pdf',collapse='',sep='')
          pdf(pdfname,paper='a4r',h=8,w=11);
          layout.m=cbind(layout.basis,layout.basis+4*(layout.basis>0));layout.m=rbind(layout.m,layout.m+8*(layout.m>0))
          layout.m[3:4,9:16]=13
          heights=c(3,1,3,1)
          plot_iteration=c(600,3000,Ni)
        }
        
        n=400
        network=make_network('SBM',n=n,clusters=10,p_within=.2,p_between=.001,rewiring=.02)
        g=network[[1]];l=network[[2]];adj=network[[3]]
        
        info_update=T
        sd_noise_information=.005
        persuasion=2
        r_min=0.1
        
        s_O=.01
        maxwell_convention=F
        
        attention_star=1
        min_attention=-.5 # to include continuous change in O as function of K
        delta_attention=0.1
        decay_attention=delta_attention/(attention_star*(n^2))
        deffuant_c=Inf
        
        # all slighty positive low attention
        information = rnorm(n, .1, 0)
        attention = runif(n, .0, .0)
        opinion=rnorm(n,0,.01)
        for( i in 1:500) opinion=stoch_cusp(n,opinion,attention+min_attention,information,s_O,maxwell_convention)
        
        # except for some negative high attention persons (activists)
        m=rep(0,n)
        m[seq(1,n,n/3)[-1]-1]=1
        information_activists=-.5
        attention_activists=1
        opinion_activists=-.5
        activism_on_at_iteration=300
        
      }
      
      if(scenario %in% c(31))
      {
        n=400
        network_sim=c('SBM','lattice')
        network=make_network(network_sim[sim_value3],n=n,clusters=10,p_within=.2,p_between=.001,rewiring=.02)
        g=network[[1]];l=network[[2]];adj=network[[3]]
        
        
        info_update=T
        Ni= 25000# 5000 # max iterations
        sd_noise_information=.0005
        persuasion=2
        r_min=0
        
        s_O=.01
        maxwell_convention=F
        
        attention_star=1
        min_attention=-.5 # to include continuous change in O as function of K
        delta_attention=0.1
        decay_attention=delta_attention/(attention_star*(n^2))
        deffuant_c=Inf
        
        delta_attention= sim_value
        r_min= sim_value2
        
        
        # all slighty positive low attention
        information = rnorm(n, .1, 0)
        attention = runif(n, .0, .0)
        opinion=rnorm(n,0,.01)
        for( i in 1:500) opinion=stoch_cusp(n,opinion,attention+min_attention,information,s_O,maxwell_convention)
        
        # except for some negative high attention persons (activists)
        m=rep(0,n)
        m[seq(1,n,n/3)[-1]-1]=1
        information_activists=-.5
        attention_activists=1
        opinion_activists=-.5
        activism_on_at_iteration=300    
        
        plot_iteration=c(1,Ni)
      }
      
      if(scenario==4) # Meat eating vegies
      {
        Ni= 30000 #150000 # max iterations
        if(pdfplot) 
        {
          set.seed(3);
          pdfname=paste('CAcusp_rev_s',scenario,'.pdf',collapse='',sep='')
          pdf(pdfname,h=8,w=11,paper='a4r'); 
          layout.m=layout.m=cbind(layout.m,layout.m+5)
          plot_iteration=c(1,Ni) 
        }
        
        n=800
        network=make_network('SBM_4',n=n,clusters=10,p_within=.2,p_between=.001)
        g=network[[1]];l=network[[2]];adj=network[[3]];gg=network[[4]]
        
        lattice=F
        info_update=T
        
        sd_noise_information=.0001
        persuasion=2
        r_min=0
        init_opinion_random = FALSE
        s_O=.01
        attention_star=1
        min_attention=-0.5 # to include continuous change in O as function of K
        delta_attention=0.01
        deffuant_c=0.2 # only interactions when |o1-o2| < .2
        p_pertubation=0.0005 # probability of inserting meat-eating vegatarians
        maxwell_convention=F
        
        information=sample(c(.1,-.4),n,T,prob=c(.8,.2))
        attention=rep(1,n)
        attention[information>0]=0.1
        opinion=rnorm(n,0,.01)
        
        information[(.5*n+1):n]=information[1:(.5*n)]
        attention[(.5*n+1):n]=attention[1:(.5*n)]
        opinion[(.5*n+1):n]=opinion[1:(.5*n)]
        
        for( i in 1:500) opinion=stoch_cusp(n,opinion,attention+min_attention,information,s_O,maxwell_convention)
        
        
        
      }
      
      if(scenario==41) # Meat eating veggies
      {
        
        n=400
        network=make_network('SBM',n=n,clusters=10,p_within=.2,p_between=.001)
        # network=make_network('lattice',n=n,clusters=10,p_within=.2,p_between=.001)
        g=network[[1]];l=network[[2]];adj=network[[3]];gg=network[[4]]
        
        lattice=F
        info_update=T
        Ni= 30000 #150000 # max iterations
        sd_noise_information=.0001
        persuasion=2
        r_min=0
        
        s_O=.01
        min_attention=-.5 # to include continuous change in O as function of K
        attention_star=1
        delta_attention=.01 # increase in attention when interacting
        deffuant_c=0.2 # only interactions when |o1-o2| < .2
        p_pertubation=0.001 # probability of inserting meat-eating vegatarians
        maxwell_convention=F
        
        p_pertubation= sim_value
        deffuant_c= sim_value2
        
        
        information=sample(c(.1,-.4),n,T,prob=c(.8,.2))
        attention=rep(1,n)
        attention[information>0]=0.1
        opinion=rnorm(n,0,.01)
        for( i in 1:500) opinion=stoch_cusp(n,opinion,attention+min_attention,information,s_O,maxwell_convention)
      }
      
      if(!pdfplot & !simulation & !PNG){layout(matrix(1:1,1,1));plot_iteration=seq(0,Ni,1000)}
      if(PNG &scenario %in% c(1,3)) plot_iteration=c(seq(1,1000,100),seq(1000,.5*Ni,500),seq(.5*Ni,Ni,1000))
      if(PNG & scenario==2) plot_iteration=c(seq(1,Ni,500))
      if(PNG &scenario %in% c(1,3,4)) plot_iteration=c(seq(1,1000,200),seq(1000,.5*Ni,1000),seq(.5*Ni,Ni,2000))
      
      
      #remember initial values
      opinion_init=opinion
      attention_init=attention
      information_init=information
      
      # make neighbor matrix for fast sampling neighbors in loop
      max_neighb=max(degree(g))
      m_neigh=matrix(NA,n,max_neighb)
      for(i in 1:n)
      {
        nb=as.numeric(neighbors(g, i, mode = c("all"))[[]])
        lengt_nb=length(nb)
        nb=c(nb,rep(NA,max_neighb-lengt_nb))
        m_neigh[i,]=nb
      }
      
      # collect data in
      if(!simulation) data=matrix(NA,Ni/200,25) else data=matrix(NA,20,25) 
      
      layout(layout.m,heights=heights)
      
      ######################################
      #############  Main loop #############
      ######################################
      
      for(iteration in 1:Ni)
      {
        if(iteration%%1000==0 & !simulation) print(iteration)
        
        # introduce activista in scenario 3 and 31 after activism_on_at_iteration interations
        if(scenario %in% c(3,31)) if(iteration == activism_on_at_iteration) # black pete opposition on
        {information[m == 1] = information_activists
        attention[m == 1] = attention_activists
        opinion[m == 1] = opinion_activists
        }
        
        #sample agent
        if(max(attention)==0) agent=0 else agent=sample(1:n,1,prob=attention)
        
        #sample partner
        if(agent!=0 & scenario!=6) 
        {
          neighb=m_neigh[agent,][!is.na(m_neigh[agent,])]
          if(length(neighb>0))  partner=sample(neighb,1) else partner=0
        } else partner=0
        if(agent!=0 & scenario==6) partner=sample(1:n,1) # in scenario 6 all agents are connected
        
        
        if(iteration>1) # exclude first iteration from updating to get a plot of initial state
        {
          
          ### information update for interacting agents
          if(partner!=0&agent!=0 ) 
          {
            I1=information[agent];A1=attention[agent];  O1=opinion[agent];
            I2=information[partner];A2=attention[partner];O2=opinion[partner]
            
            if( abs(O1-O2) < deffuant_c) # bounded confidence check
            {
              information[agent]=information_update(I1,I2,A1,A2,persuasion,r_min,info_update)
              information[partner]=information_update(I2,I1,A2,A1,persuasion,r_min,info_update)
            }
            
            ### attention update for interacting agents
            if( abs(O1-O2) < deffuant_c )
            {
              attention[agent]=attention[agent]+delta_attention*(2*attention_star-attention[agent]) 
              attention[partner]=attention[partner]+delta_attention*(2*attention_star-attention[partner])
            }
          }
          
          ### attention decay for all agents
          information=information+rnorm(n,0,sd_noise_information)
          attention=attention-2*delta_attention*attention/n   # correction 2 times if interaction
          
          
          if(scenario==2) if(iteration>(1/3)*Ni) {information=.999*information;sd_noise_information=0}  # schrinking I in scenario 1
          if(scenario==2) if(iteration>(2/3)*Ni) {attention=.999*attention;delta_attention=0} # schrinking A in scenario 1
          
          ###  behavior update
          opinion=stoch_cusp(n,opinion,attention+min_attention,information,s_O,maxwell_convention)
          
          if(scenario==4) # introducing meat eating vegatations in scenario 4
          {
            perturb_vector=sample(c(F,T),n,T,prob=c(1-p_pertubation,p_pertubation))
            opinion[opinion<0 & perturb_vector & rep(c(T,F),each=.5*n)] = -1*opinion[opinion<0 & perturb_vector & rep(c(T,F),each=.5*n)]
          }
          
          
          if(scenario==41) # introducing meat eating vegatations in scenario 4
          {
            perturb_vector=sample(c(F,T),n,T,prob=c(1-p_pertubation,p_pertubation))
            opinion[opinion<0 & perturb_vector] = -1*opinion[opinion<0 & perturb_vector]
          }
          
        }
        
        # reports
        if(partner!=0&agent!=0 & iteration %% 200==0 & !simulation)  
        {
          CD=27*information^2-4*(attention+min_attention)^3
          CD_p=sum(CD<0)/n
          freq_pos_opinion=sum(opinion>0)/n
          assortativity_g=assortativity(g,opinion)
          ambivalence_r=cor(opinion,information)
          dip=dip(opinion)
          
          data[iteration/200,]=c(agent,partner,O1,O2,I1,I2,A1,A2,opinion[agent],opinion[partner],information[agent],
                                 information[partner],attention[agent],attention[partner],
                                 CD_p,ambivalence_r,freq_pos_opinion,
                                 assortativity_g,dip,mean(opinion),sd(opinion),mean(information),
                                 sd(information),mean(attention),sd(attention))
        }
        
        if(partner!=0&agent!=0 & iteration > (Ni-21) & simulation)  
        {
          CD=27*information^2-4*(attention+min_attention)^3
          CD_p=sum(CD<0)/n
          freq_pos_opinion=sum(opinion>0)/n
          assortativity_g=assortativity(g,opinion)
          ambivalence_r=cor(opinion,information)
          dip=dip(opinion)
          
          data[iteration-(Ni-20),]=c(agent,partner,O1,O2,I1,I2,A1,A2,opinion[agent],opinion[partner],information[agent],
                                     information[partner],attention[agent],attention[partner],
                                     CD_p,ambivalence_r,freq_pos_opinion,
                                     assortativity_g,dip,mean(opinion),sd(opinion),mean(information),
                                     sd(information),mean(attention),sd(attention))
        }
        
        
        ###### Plotting
        if(plotting)
        {
          if(PNG) {
            png(paste0("pngplots_",scenario,"/",scenario,"_",100000+iteration,".png"),width=8,height=5,units="in",res=200)
            layout(layout.m,heights=heights)
          }
          
          if(iteration %in% plot_iteration)
          {
            min.o=-1.5;max.o=1.5;min.i=-1;max.i=1;min.a=0;max.a=2
            opin=sapply(-opinion, function(y) min(max(y,min.o),max.o))
            inform=sapply(-information, function(y) min(max(y,min.i),max.i))
            atten=attention;atten[atten>2]=2
            CD=27*information^2-4*(attention+min_attention)^3
            assortativity_g=assortativity(g,opinion)
            
            shape=rep("circle",n)
            if(pdfplot)c.title=1.1  else c.title=1.5
            title=sub_H=sub_A=""
            if(scenario==2 )
            {
              h_d=hartigan_d(opinion)
              if (iteration >= 1) title = "Initial state"
              if (iteration >= (1/3)*Ni) title = "All A = 1, I varies"
              if (iteration >= (2/3)*Ni) title= "A = 1, I = 0"
              if (iteration == Ni) title = "A = 0, I = 0"
              sub_H=paste('\nHartigan D = ',round(h_d[[1]],2),h_d[[2]],sep="",col="")
              sub_A=paste("\nAssortativity =",round(assortativity_g,2))
              title=paste(title,sub_H,sub_A)
              sub=""
            }
            
            if(scenario==3 )
            {
              h_d=hartigan_d(opinion)
              title=paste('p(O>0) =',round(sum(opinion>0)/n,2),'\nsd(O) =',round(sd(opinion),2))
              shape[m==1]='star'
              sub_H=paste('\nHartigan D = ',round(h_d[[1]],2),h_d[[2]],sep="",col="")
              sub_A=paste("\nAssortativity =",round(assortativity_g,2))
              title=paste(title,sub_H,sub_A)
            }
            
            if(scenario %in% c(2,3)) 
            {
              plot.graph(adj,l,opin,inform,atten,CD,title,shape,c.title)
              plot.histo(-opin,min.o,max.o,'O')
              plot.histo(-inform,min.i,max.i,'I')
              plot.histo(atten, min.a,max.a,'A')
            }
            
            if(scenario==4){
              
              h_d=hartigan_d(opinion[(1+(n/2)):n])
              title=paste('p(O>0) =',round(sum(opinion[(1+(n/2)):n]>0)/(.5*n),2),'\nsd(O) =',round(sd(opinion[(1+(n/2)):n]),2))
              sub_H=paste('\nHartigan D = ',round(h_d[[1]],2),h_d[[2]],sep="",col="")
              sub_A=paste("\nAssortativity =",round(assortativity( gg,opinion[(1+(n/2)):n]),2))
              sub1=paste(title,sub_H,sub_A)
              
              h_d=hartigan_d(opinion[1:(n/2)])
              sub=paste('p(O>0) =',round(sum(opinion[1:(n/2)]>0)/(.5*n),2),'\nsd(O) =',round(sd(opinion[1:(n/2)]),2))
              sub_H=paste('\nHartigan D = ',round(h_d[[1]],2),h_d[[2]],sep="",col="")
              sub_A=paste("\nAssortativity =",round(assortativity(gg,opinion[1:(n/2)]),2))
              sub=paste(sub,sub_H,sub_A)
              
              plot.graph(adj,l,opin,inform,atten,CD,title="",shape,c.title)
              par(mar=c(0,0,0,0))
              plot(c(0,1),c(0,1),type='n',axes=FALSE,ann=FALSE)
              text(0,.3,sub1,cex=.9,pos=4,family='mono')
              plot.histo(-opin[(1+(n/2)):n],min.o,max.o,'O')
              par(mar=c(0,0,0,0))
              plot(c(0,1),c(0,1),type='n',axes=FALSE,ann=FALSE)
              text(0,.3,sub,cex=.9,pos=4,family='mono')
              plot.histo(-opin[1:(n/2)],min.o,max.o,'O')
            }  
            
            if(scenario==1 )
            {
              #image(opinioni,axes=F,col = heat.colors(256),breaks=seq(-2,2,le=257),xlab='I', ylab='K(init)')
              opin_m=matrix(opin,n1,n2)
              image(opin_m,axes=F,col = diverging_hsv(256,s=1.6,v=1),breaks=seq(-1.5,1.5,le=257),xlab='A', ylab='I(init)')
              axis(1,seq(0,1,length=6),seq(i_init_min,i_init_max,by=.1))
              axis(2,seq(0,1,by=.25),seq(-1*k_init,k_init,length=5))
              scale_factor= (i_init_max-i_init_min)/(2*k_init)
              bf2=function(b) .5+scale_factor*(2 * (b+min_attention+i_init_min)^(3/2))/(3 *3^(1/2))
              plot(bf2,0,1,ylab='b',xlab='t',add=T)
              bf2=function(b) .5-1*scale_factor*(2 * (b+min_attention+i_init_min)^(3/2))/(3 *3^(1/2))
              plot(bf2,0,1,ylab='b',add=T)
            }
            
          }      
          if(PNG) dev.off()
        }
      }
      ####################
      ########end#########
      ####################
      
      # cusp fig plot in scenario 3
      if(scenario==3 & pdfplot) 
      {cuspfig <- readPNG('figure7.png')
      plot(c(0,1),c(0,1), type='n', main="", xlab="", ylab="",axes=F)
      par(mar=c(0,0,0,0))
      rasterImage(cuspfig, .05,-.05, .9,1.05)
      }
      
      if(pdfplot) dev.off()
      
      datad=as.data.frame(data)
      names(datad)=c('agent','partner','o1o','o2o','i1o','i2o','a1o','a2o','o1n','o2n','i1n','i2n','a1n','a2n','CD',
                     'ambivalence_r','freq_pos_opinion','assortativity','dip','mean_o','sd_o',
                     'mean_i','sd_i','mean_a','sd_a')
      
      datasim[sim_i,]=c(sim_value,sim_value2,sim_value3,mean(datad$freq_pos_opinion),mean(datad$sd_o),
                        mean(datad$dip),mean(datad$assortativity))
    }
  }
} # end of simulation loops

if(scenario==31)
{
  
  mean_sim=as.matrix(aggregate(datasim[,4:7],list(datasim[,1],datasim[,2],datasim[,3]),mean,na.rm=T))
  colnames(mean_sim)=c(sim_var,sim_var2,sim_var3,c('p(O>0)','SD(O)','Hartigan D','Assortativity'))
  mean_sim2=as_tibble(mean_sim)
  mean_sim2$r_min=as.factor(mean_sim2$r_min)
  mean_sim2$network[mean_sim2$network==1]='Stoch. block model'
  mean_sim2$network[mean_sim2$network==2]='Lattice network'
  mean_sim2$network=factor(mean_sim2$network, levels=c('Stoch. block model','Lattice network'))
  
  mean_sim3 = mean_sim2 %>%
    gather('p(O>0)','Hartigan D', key = measure, value = value)
  
  mean_sim3$measure=factor(mean_sim3$measure, levels=c('p(O>0)','Hartigan D'))
  
  
  #save(datasim,file='datasim_sc31')
  pdf('CAcusp_rev_s31.pdf')
  
  ggplot(data = mean_sim3) + 
    geom_point(mapping = aes(x = d_A, y = value,color = r_min),size=1.2)+
    geom_line(mapping = aes(x = d_A, y = value,color = r_min),size=.5)+
    facet_wrap(network~measure,scales="free", dir = "h")+
    scale_color_grey() + theme_classic()+ xlab(expression('d'['A']))+ylab('')+
    guides(color = guide_legend(reverse = TRUE))+
    theme(strip.background = element_rect(fill=NA,colour = NA))+
    theme(strip.text.x = element_text(size = 12))+
    theme(strip.text.x = element_text(margin = margin(1, 0, 1, 0)))+
    labs(color = expression('r'['min']))
  
  
  ggplot(data = mean_sim3) + 
    geom_point(mapping = aes(x = d_A, y = value,color = r_min),size=1.2)+
    geom_line(mapping = aes(x = d_A, y = value,color = r_min),size=.5)+
    facet_wrap(network~measure,scales="free", dir = "h")+
    scale_colour_hue() + theme_classic()+ xlab(expression('d'['A']))+ylab('')+
    guides(color = guide_legend(reverse = TRUE))+
    theme(strip.background = element_rect(fill=NA,colour = NA))+
    theme(strip.text.x = element_text(size = 12))+
    theme(strip.text.x = element_text(margin = margin(1, 0, 1, 0)))+
    labs(color = expression('r'['min']))
  
  
  dev.off()
  
}

if(scenario==41)
{
  
  mean_sim=as.matrix(aggregate(datasim[,4:7],list(datasim[,1],datasim[,2],datasim[,3]),mean,na.rm=T))
  colnames(mean_sim)=c('p_perturbation',sim_var2,sim_var3,c('p(O>0)','SD(O)','Hartigan D','Assortativity'))
  mean_sim2=as_tibble(mean_sim)
  mean_sim2$t_d=as.factor(mean_sim2$t_d)
  mean_sim2$network[mean_sim2$network==1]='Stoch. block model'
  mean_sim2$network[mean_sim2$network==2]='Lattice network'
  mean_sim2$network=factor(mean_sim2$network, levels=c('Stoch. block model','Lattice network'))
  
  mean_sim3 = mean_sim2 %>%
    gather('p(O>0)','Hartigan D', key = measure, value = value)
  mean_sim3$measure=factor(mean_sim3$measure, levels=c('p(O>0)','Hartigan D'))
  
  
  #  save(datasim,file='datasim_sc41')
  pdf('CAcusp_rev_s41.pdf')
  
  
  scaleFUN <- function(x) sprintf("%.4f", x)
  
  ggplot(data = mean_sim3) + 
    geom_point(mapping = aes(x = p_perturbation, y = value,color = t_d),size=1.2)+
    geom_line(mapping = aes(x = p_perturbation, y = value,color = t_d),size=.5)+
    facet_wrap(network~measure,scales="free", dir = "h")+
    scale_color_grey() + theme_classic()+ xlab('p(pertubation)')+ylab('')+
    guides(color = guide_legend(reverse = TRUE))+
    theme(strip.background = element_rect(fill=NA,colour = NA))+
    theme(strip.text.x = element_text(size = 12))+
    theme(strip.text.x = element_text(margin = margin(2, 2, 2, 0)))+
    labs(color = expression('t'['O']))
  
  ggplot(data = mean_sim3) + 
    geom_point(mapping = aes(x = p_perturbation, y = value,color = t_d),size=1.2)+
    geom_line(mapping = aes(x = p_perturbation, y = value,color = t_d),size=.5)+
    facet_wrap(network~measure,scales="free", dir = "h")+
    scale_colour_hue() + theme_classic()+ xlab('p(pertubation)')+ylab('')+
    guides(color = guide_legend(reverse = TRUE))+
    theme(strip.background = element_rect(fill=NA,colour = NA))+
    theme(strip.text.x = element_text(size = 14))+
    theme(strip.text.x = element_text(margin = margin(1, 0, 1, 0)))+
    labs(color = expression('t'['O']))
  
  
  dev.off()
  
}

layout(1:3)
hist(opinion,20,col='grey');hist(information,20,col='grey');hist(attention,20,col='grey')

layout(1)
matplot(datad[,c('CD','freq_pos_opinion','assortativity','ambivalence_r')],type='l',bty='n',xlab='x 100',ylab="value",axes=F,col=1:4,lty=1:4)
legend('topright',col=1:7,lty=1:7,legend=c('in bifurcation set','ambivalence','freq_pos_opinion','assortativity','ambivalence_r'),bty='n')
axis(2)
axis(1)

if(PNG)
{
  library(dplyr)
  library(purrr) 
  library(magick)
  
  list.files(path=paste0("pngplots_",scenario,"/"), pattern = '*.png', full.names = TRUE) %>% 
    image_read() %>% # reads each path file
    image_join() %>% # joins image
    image_animate(fps=2,loop=0) %>% # animates, can opt for number of loops
    image_write(paste0("Anim_",scenario,".gif")) # write to current dir
  
}

# install.packages('plotly')
# library(plotly)
# 
# plot_ly(x = ~as.vector(information), y = ~as.vector(attention), z = ~as.vector(opinion),
#              marker = list(size=1+attention*4,color = ~as.vector(-opinion), colorscale = c('#FF0000FF', '#FFFFFDFF'), showscale = TRUE)) %>%
#   add_markers() %>%
#   layout(scene = list(xaxis = list(title = 'Information'),
#                       yaxis = list(title = 'Attention'),
#                       zaxis = list(title = 'Opinion')),
#          annotations = list(
#            x = 1.13,
#            y = 1.05,
#            text = 'Polarization',
#            xref = 'paper',
#            yref = 'paper',
#            showarrow = FALSE
#          ))
# 
# library(cusp)
# ?cusp
# ?cusp3d.surface
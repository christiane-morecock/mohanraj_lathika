# Version dated 18/06-2014.
# This version corrects a programmeing error in line 121.

Normfinder=function(filename,Groups=TRUE,ctVal=TRUE,pStabLim=0.25){
  #
  # If Groups is TRUE the last row contains the group identifier,
  # and the last row must be started by a name for that row.
  # No spaces are allowed in the gene names, sample names and group identifier.
  #
  dat0=read.table(filename,header=TRUE,row.names=1,colClasses="character")
  # dat0 = filename
  #
  ntotal=dim(dat0)[2] # number of samples
  k0=dim(dat0)[1] # number of rows
  #
  if (Groups){
    ngenes=k0-1 # number of genes
    genenames=rownames(dat0)[-k0]
    grId=dat0[k0,]
    dat0=dat0[-k0,]
  } else {
    ngenes=k0 # number of genes
    genenames=rownames(dat0)
    grId=rep(1,ntotal)
  }
  #
  dat=matrix(as.numeric(unlist(dat0)),ngenes,ntotal) # matrix instead of list
  #
  if (!ctVal){dat=log2(dat)} # transform to log2 values
  #
  samplenames=colnames(dat0)
  grId=factor(unlist(grId))  # group identifier
  groupnames=levels(grId)  # group names
  ngr=length(levels(grId)) # number of groups
  # Number of samples in each group:
  nsamples=rep(0,ngr)
  for (group in 1:ngr){nsamples[group]=sum(grId==groupnames[group])}
  #
  #
  MakeStab=function(da){
    ngenes=dim(da)[1]
    # Sample averages
    sampleavg=apply(da,2,mean)
    # Gene averages within group
    genegroupavg=matrix(0,ngenes,ngr)
    for (group in 1:ngr){
      genegroupavg[,group]=apply(da[,grId==groupnames[group]],1,mean)}
    # Group averages
    groupavg=rep(0,ngr)
    for (group in 1:ngr){groupavg[group]=mean(da[,grId==groupnames[group]])}
    
    # Variances 
    GGvar=matrix(0,ngenes,ngr)
    for (group in 1:ngr){
      grset=(grId==groupnames[group])
      a=rep(0,ngenes)
      for (gene in 1:ngenes){
        a[gene]=sum((da[gene,grset]-genegroupavg[gene,group]-
                       sampleavg[grset]+groupavg[group])^2)/(nsamples[group]-1)
      }
      GGvar[,group]=(a-sum(a)/(ngenes*ngenes-ngenes))/(1-2/ngenes)
    }
    #
    # Change possible negative values
    genegroupMinvar=matrix(0,ngenes,ngr)
    for (group in 1:ngr){
      grset=(grId==groupnames[group])
      z=da[,grset]
      for (gene in 1:ngenes){
        varpair=rep(0,ngenes)
        for (gene1 in 1:ngenes){varpair[gene1]=var(z[gene,]-z[gene1,])}
        genegroupMinvar[gene,group]=min(varpair[-gene])/4
      }
    }
    #
    # Final variances
    GGvar=ifelse(GGvar<0,genegroupMinvar,GGvar)
    #
    # Old stability measure for each gene is calculated:
    #
    dif=genegroupavg
    difgeneavg=apply(dif,1,mean)
    difgroupavg=apply(dif,2,mean)
    difavg=mean(dif)
    for (gene in 1:ngenes){
      for (group in 1:ngr){
        dif[gene,group]=dif[gene,group]-difgeneavg[gene]-difgroupavg[group]+difavg
      }
    }
    #
    nsampMatrix=matrix(rep(nsamples,ngenes),ngenes,ngr,byrow=T)
    vardif=GGvar/nsampMatrix
    gamma=sum(dif*dif)/((ngr-1)*(ngenes-1))-sum(vardif)/(ngenes*ngr)
    gamma=ifelse(gamma<0,0,gamma)
    #
    difnew=dif*gamma/(gamma+vardif)
    varnew=vardif+gamma*vardif/(gamma+vardif)
    Ostab0=abs(difnew)+sqrt(varnew)
    Ostab=apply(Ostab0,1,mean)
    #
    # Measure of group differences:
    mud=rep(0,ngenes)
    for (gene in 1:ngenes){
      mud[gene]=2*max(abs(dif[gene,]))
    }
    # Common variance:
    genevar=rep(0,ngenes)
    for (gene in 1:ngenes){
      genevar[gene]=sum((nsamples-1)*GGvar[gene,])/(sum(nsamples)-ngr)
    }
    Gsd=sqrt(genevar)
    #
    # Return results:
    #
    return(cbind(mud,Gsd,Ostab,rep(gamma,ngenes),GGvar,dif))
  }    # End of function MakeStab
  #
  #
  MakeComb2=function(g1,g2,res){
    gam=res[1,4]
    d1=res[g1,(4+ngr+1):(4+ngr+ngr)]; d2=res[g2,(4+ngr+1):(4+ngr+ngr)]
    s1=res[g1,(4+1):(4+ngr)]^2; s2=res[g2,(4+1):(4+ngr)]^2
    rho=abs(gam*d1/(gam+s1/nsamples)+gam*d2/(gam+s2/nsamples))*
      sqrt(ngenes/(ngenes-2))/2
    rho=rho+sqrt(s1/nsamples+gam*s1/(nsamples*gam+s1)+
                   s2/nsamples+gam*s2/(nsamples*gam+s2))/2
    return(sum(rho)/2)
  }
  #
  #
  MakeStabOne=function(da){
    ngenes=dim(da)[1]
    # Sample averages
    sampleavg=apply(da,2,mean)
    # Gene averages
    geneavg=apply(da,1,mean)
    totalavg=mean(da)
    #
    # Variances 
    genevar0=rep(0,ngenes)
    for (gene in 1:ngenes){
      genevar0[gene]=
        sum((dat[gene,]-geneavg[gene]-sampleavg+totalavg)^2)/
        ((ntotal-1)*(1-2/ngenes))
    }
    genevar=genevar0-sum(genevar0)/(ngenes*ngenes-ngenes)
    #
    # Change possible negative values
    geneMinvar=rep(0,ngenes)
    z=da
    for (gene in 1:ngenes){
      varpair=rep(0,ngenes)
      for (gene1 in 1:ngenes){varpair[gene1]=var(z[gene,]-z[gene1,])}
      geneMinvar[gene]=min(varpair[-gene])/4
    }
    # Final variances
    genevar=ifelse(genevar<0,geneMinvar,genevar)
    #
    return(genevar)
  }
  #     End of function MakeStabOne
  #
  #################################################
  #
  # Main part
  #
  if (ngr>1){   # More than one group.
    #
    res=MakeStab(dat)
    #
    gcand=c(1:ngenes)[res[,3]<pStabLim]
    ncand=length(gcand)
    if (ncand<4){
      if (ngenes>3){
        li=sort(res[,3])[4]
        gcand=c(1:ngenes)[res[,3]<=li]
        ncand=length(gcand)
      } else {
        gcand=c(1:ngenes)
        ncand=length(gcand)
      }
    }
    #
    vv2=c()
    #
    for (g1 in 1:(ncand-1)){
      for (g2 in (g1+1):ncand){
        qmeas=MakeComb2(gcand[g1],gcand[g2],res)
        vv2=rbind(vv2,c(gcand[g1],gcand[g2],qmeas))
      }}
    #
    ord=order(res[,3])
    FinalRes=list(Ordered=
                    data.frame("GroupDif"=round(res[ord,1],2),"GroupSD"=round(res[ord,2],2),
                               "Stability"=round(res[ord,3],2),row.names=genenames[ord]),
                  UnOrdered=
                    data.frame("GroupDif"=round(res[,1],2),"GroupSD"=round(res[,2],2),
                               "Stability"=round(res[,3],2),
                               "IGroupSD"=round(sqrt(res[,(4+1):(4+ngr)]),2),
                               "IGroupDif"=round(res[,(4+ngr+1):(4+ngr+ngr)],2),
                               row.names=genenames),
                  PairOfGenes=
                    data.frame("Gene1"=genenames[vv2[,1]],"Gene2"=genenames[vv2[,2]],
                               "Stability"=round(vv2[,3],2)))
    #
    return(FinalRes)
    #
  } else {    # End of more than one group: next is for one group only.
    #
    #
    sigma=sqrt(MakeStabOne(dat))
    #
    siglim=(min(sigma)+0.1)
    gcand=c(1:ngenes)[sigma<siglim]
    ncand=length(gcand)
    #
    if ((ncand>=2)&(ngenes>3)){
      #
      vv2=c()
      #
      for (g1 in 1:(ncand-1)){
        for (g2 in (g1+1):ncand){
          dat1=rbind(dat[-c(gcand[g1],gcand[g2]),],
                     apply(dat[c(gcand[g1],gcand[g2]),],2,mean))
          qmeas=sqrt(MakeStabOne(dat1))
          vv2=rbind(vv2,c(gcand[g1],gcand[g2],qmeas[ngenes-1]))
        }}
      ord=order(sigma)
      FinalRes=list(Ordered=
                      data.frame("GroupSD"=round(sigma[ord],2),row.names=genenames[ord]),
                    PairOfGenes=
                      data.frame("Gene1"=genenames[vv2[,1]],"Gene2"=genenames[vv2[,2]],
                                 "GroupSD"=round(vv2[,3],2)))
    } else { # No combined genes to consider
      ord=order(sigma)
      FinalRes=list(Ordered=
                      data.frame("GroupSD"=round(sigma[ord],2),row.names=genenames[ord]))
    } # End ncand<2 or ngenes<=3
    #
    return(FinalRes)
    #
  }  # End one group only
  #
} # End of main function

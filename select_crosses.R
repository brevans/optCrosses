setwd("E:/collaborations/EEB/201312/")
x<-read.table("crosses_HiRes.txt",sep='\t',header=T)
y<-x[,2:(ncol(x)-3)]
y<-as.matrix(y)
recode<-lapply(1:nrow(y),function(ii){
  gt<-numeric(ncol(y))
  gt[y[ii,]=="BB"|y[ii,]=="AA"]<-0
  gt[y[ii,]=="AB"]<-1
  gt[y[ii,]=="NoCall"]<-(-9)
  gt
})
recode<-matrix(unlist(recode),nrow=nrow(y),byrow=T)
crosses<-c(1:(ncol(y)/2))
count.info<-function(crossi){
  gtmat<-recode[,c(2*crossi-1,2*crossi)]
  info.col<-numeric(nrow(gtmat))
  info.col[gtmat[,1]+gtmat[,2]==1]<-1
  info.col
}
cross.info<-lapply(crosses,FUN="count.info")
cross.info.mat<-matrix(unlist(cross.info),ncol=length(crosses),byrow=F)
system.time(
for(i1 in 1:23)
  for(i2 in (i1+1):24)
    for(i3 in (i2+1):25)
      for(i4 in (i3+1):26)
        for(i5 in (i4+1):27)
         cat(i1,'\t',i2,'\t',i3,'\t',i4,'\t',i5,'\t',
          sum(cross.info.mat[,i1]+cross.info.mat[,i2]+cross.info.mat[,i3]+cross.info.mat[,i4]+cross.info.mat[,i5]>=1),'\n',
          file="cross_info_5.txt",append=T)
)
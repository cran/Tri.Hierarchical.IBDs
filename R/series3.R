#' Tri-Hierarchical IBDs using Rectangular Association Scheme
#'
#' @param m Any integer >=3
#' @param n Any integer >=3
#' @param D1 Bi-Hierarchical IBD by ignoring blocks
#' @param D2 Bi-Hierarchical IBD by ignoring sub-blocks
#' @param D3 Bi-Hierarchical IBD by ignoring sub-sub blocks
#' @param D4 IBD at block level
#' @param D5 IBD at sub block level
#' @param D6 IBD at sub-sub block level
#' @param Randomization Randomization of layout of the designs if needed enter TRUE; by default it is FALSE.
#'@description
#'This function provides the Tri-Hierarchical IBDs based on Rectangular association scheme. Here,  v= m*n, v should be composite number and (m,n)>=3. We find balanced incomplete block designs (BIBD) at block level and rectangular PBIB designs at sub-block level as well as sub-sub block level. Information matrix pertaining to the estimation of treatments effects, canonical efficiency factor in comparison to an orthogonal design and six component designs are provided.
#'
#' @return It gives Tri-HIB design and six component designs with canonical efficiency factor in comparison to an orthogonal design.
#' @export
#'@note Numbers in the outer most parentheses represents as block elements, second level parentheses as sub block elements and inner most parentheses as sub-sub block elements.
#'
#'@references
#'Preece, D.A. (1967) <https://doi.org/10.1093/biomet/54.3-4.479>."Nested balanced incomplete block designs".
#' @examples
#' library(Tri.Hierarchical.IBDs)
#' Series3(4,3,D1=TRUE,D2=TRUE, D3=TRUE, D4=TRUE,D5=FALSE,D6=TRUE,Randomization=TRUE)
Series3=function(m,n,D1=FALSE,D2=FALSE,D3=FALSE, D4=FALSE,D5=FALSE,D6=FALSE,Randomization=FALSE){
if(m>=3 && n>=3){
  if(m==3 && n==3){
    message("Please enter m and n such that v = m*n and v>=12")
  }else{
  v=m*n
################
IBD_Cannonical=function(Design){
  ysd=Design
  length_ysd=length(ysd)-length(ysd[ysd==0])
  #######obsn vs trt
  delprime=NULL
  for(i in 1:nrow(ysd)){
    for(j in 1:ncol(ysd)){
      if(ysd[i,j]>0){
        create_vec=matrix(0,nrow=1,ncol=max(ysd))
        ele=ysd[i,j]
        create_vec[,ele]<-1
        delprime=rbind(delprime,create_vec)
      }else{
        j=j+1
      }
    }
  }
  length(ysd[ysd>0])
  #######
  rep_mat= t(delprime)%*%delprime
  ######################################### observation vs block matrix
  ################################################
  D1_mat_prime=NULL
  for(i in 1:nrow(ysd)){
    rowsize=length(ysd[i,][ysd[i,]>0])
    zeromatrix=matrix(0,nrow=rowsize,ncol=nrow(ysd))
    zeromatrix[,i]<-1
    D1_mat_prime=rbind(D1_mat_prime,zeromatrix)
  }
  K_matrix=t(D1_mat_prime)%*%D1_mat_prime
  N_matrix=t(delprime)%*%D1_mat_prime
  #library(MASS)
  c_mat=rep_mat-N_matrix%*%solve(K_matrix)%*%t(N_matrix)
  e1=eigen(c_mat)$values
  e1=e1[e1>0.000000001]
  e2=e1/rep_mat[1,1]
  e3=1/e2
  cefficiency=length(e3)/sum(e3)
  cefficiency=round(cefficiency,digits = 4)
  l1=list("Cannonical Efficiency Factor"=cefficiency)
  print(l1)
}


######################
recmat=matrix(1:v,nrow=m,ncol=n,byrow=T)
final=NULL
for(i in recmat){
  vec=which(recmat==i,arr.ind = T)
  row=vec[1]
  col=vec[2]
  AS1=setdiff(recmat[row,],i)
  AS2=setdiff(recmat[,col],i)
  AS3=setdiff(recmat,c(i,AS1,AS2))
  final=rbind(final,c(AS1,AS2,AS3))
}
#write.table(final,"clipboard",sep="\t")
n1=length(AS1)
n2=length(AS2)
n3=length(AS3)
###########################
rand=c(sample(1:nrow(final)))
if(Randomization==T){
  final=final[c(rand),]
}
##########################
final1=cbind("{","[","(",final[,1:n1],")","(",final[,(n1+1):(n1+n2)],")","]","[","(",final[,(1+n2+n1):ncol(final)],")","]","}")
blk=c()
for(i in 1:nrow(final)){
  blk=c(blk,paste0("Block-",i,collapse = ""))
}
row.names(final1)<-blk
message("Tri-Hierarchical IBD")
cat("\n")
prmatrix(final1, rowlab = , collab = rep("", ncol(final1)),quote = F)
#######################################################
v=m*n
r=v-1
cat("\n")
list1=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),
        paste("Number of  blocks(b1) =",v),
        paste("Number of sub blocks(b2) =",2*v), paste("Number of sub-sub blocks(b3) =",3*v),
        paste("Block size(k1) =",v-1),
        paste("1st sub block size(k21) =",(n1+n2)),
        paste("2nd sub block size(k22) =",n3),
        paste("1st sub-sub block size(k31) =",n1),
        paste("2nd sub-sub block size(k32) =",n2),
        paste("2nd sub-sub block size(k33) =",n3),
        paste("Number of times each pair occur together in blocks (lambda_1) =",(m*n-2)),
        paste("Number of times pair of 1st associates occur together in sub blocks (lambda_21) =",m*(n-2)),
        paste("Number of times pair of 2nd associates occur together in sub blocks (lambda_22) =",n*(m-2)),
        paste("Number of times pair of 3rd associates occur together in sub blocks (lambda_23) =",m*(n-2)),
        paste("Number of times pair of 1st associates occur together in sub-sub blocks (lambda_31) =",m*(n-2)),
        paste("Number of times pair of 2nd associates occur together in sub-sub blocks (lambda_32) =",n*(m-2)),
        paste("Number of times pair of 3rd associates occur together in sub-sub blocks (lambda_33) =",(n-2)*(m-2)))
print(list1,quote=F)
############eff mat
eff_mat1=final[,1:n1]
eff_mat2=final[,(n1+1):(n1+n2)]
eff_mat3=final[,(n1+n2+1):(n1+n2+n3)]
gap1=ncol(eff_mat3)-ncol(eff_mat1)
gap2=ncol(eff_mat3)-ncol(eff_mat2)
eff_mat1=cbind(eff_mat1,matrix(0,nrow=nrow(eff_mat1),ncol=gap1))
eff_mat2=cbind(eff_mat2,matrix(0,nrow=nrow(eff_mat2),ncol=gap2))
new_eff=rbind(eff_mat1,eff_mat2,eff_mat3)
IBD_Cannonical(new_eff)
cat("\n")
if(D1==T){
  message("Bi-Hierarchical IBD by ignoring blocks")
  cat("\n")
  remove_vec=c(1,ncol(final1))
  prmatrix(final1[,-c(remove_vec)], rowlab = , collab = rep("", ncol(final1[,-c(remove_vec)])),quote = F)
  cat("\n")
  ###eff mat
  eff_mat1=final[,1:n1]
  eff_mat2=final[,(1+n1):(n1+n2)]
  eff_mat3=final[,(n1+n2+1):(ncol(final))]
  gap1=ncol(eff_mat3)-ncol(eff_mat1)
  gap2=ncol(eff_mat3)-ncol(eff_mat2)
  eff_mat1=cbind(eff_mat1,matrix(0,nrow=nrow(eff_mat1),ncol=gap1))
  eff_mat2=cbind(eff_mat2,matrix(0,nrow=nrow(eff_mat2),ncol=gap2))
  new_eff=rbind(eff_mat1,eff_mat2,eff_mat3)
  ###Parameters need to add
  list1=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),
          paste("Number of sub blocks(b2) =",2*v), paste("Number of sub-sub blocks(b3) =",3*v),
          paste("1st sub block size(k21) =",(n1+n2)),paste("2nd sub block size(k22) =",n3),
          paste("1st sub-sub block size(k31) =",n1),paste("2nd sub-sub block size(k32) =",n2),paste("2nd sub-sub block size(k33) =",n3),
          paste("Number of times pair of 1st associates occur together in sub blocks (lambda_21) =",m*(n-2)),
          paste("Number of times pair of 2nd associates occur together in sub blocks (lambda_22) =",n*(m-2)),
          paste("Number of times pair of 3rd associates occur together in sub blocks (lambda_23) =",m*(n-2)),
          paste("Number of times pair of 1st associates occur together in sub-sub blocks (lambda_31) =",m*(n-2)),
          paste("Number of times pair of 2nd associates occur together in sub-sub blocks (lambda_32) =",n*(m-2)),
          paste("Number of times pair of 3rd associates occur together in sub-sub blocks (lambda_33) =",(n-2)*(m-2))
          )
  print(list1,quote=F)
  IBD_Cannonical(new_eff)
  cat("\n")
}
if(D2==T){
  message("Bi-Hierarchical IBD by ignoring sub blocks")
  cat("\n")
  remove_vec=c(2,(3+n1+2+n2+1+1),(3+n1+2+n2+1+2),ncol(final1)-1)
  prmatrix(final1[,-c(remove_vec)], rowlab = , collab = rep("", ncol(final1[,-c(remove_vec)])),quote = F)
  cat("\n")
  ###eff mat
  eff_mat1=final[,1:n1]
  eff_mat2=final[,(1+n1):(n1+n2)]
  eff_mat3=final[,(n1+n2+1):(ncol(final))]
  gap1=ncol(eff_mat3)-ncol(eff_mat1)
  gap2=ncol(eff_mat3)-ncol(eff_mat2)
  eff_mat1=cbind(eff_mat1,matrix(0,nrow=nrow(eff_mat1),ncol=gap1))
  eff_mat2=cbind(eff_mat2,matrix(0,nrow=nrow(eff_mat2),ncol=gap2))
  new_eff=rbind(eff_mat1,eff_mat2,eff_mat3)
  ###Parameters need to add
  list1=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),
          paste("Number of  blocks(b1) =",v),
          paste("Number of sub-sub blocks(b3) =",3*v),
          paste("Block size(k1) =",v-1),paste("1st sub block size(k21) =",(n1+n2)),paste("2nd sub block size(k22) =",n3),
          paste("1st sub-sub block size(k31) =",n1),paste("2nd sub-sub block size(k32) =",n2),paste("2nd sub-sub block size(k33) =",n3),
          paste("Number of times each pair occur together in blocks (lambda_1) =",(m*n-2)),
                paste("Number of times pair of 1st associates occur together in sub-sub blocks (lambda_31) =",m*(n-2)),
                paste("Number of times pair of 2nd associates occur together in sub-sub blocks (lambda_32) =",n*(m-2)),
                paste("Number of times pair of 3rd associates occur together in sub-sub blocks (lambda_33) =",(n-2)*(m-2))
          )
  print(list1,quote=F)
  IBD_Cannonical(new_eff)
  cat("\n")

}
if(D3==T){
  message("Bi-Hierarchical IBD by ignoring sub-sub blocks")

  cat("\n")
  remove_vec=c(3,(3+n1+1),(3+n1+2),(3+n1+2+n2+1),(3+n1+2+n2+1+2+1),ncol(final1)-2)
  prmatrix(final1[,-c(remove_vec)], rowlab = , collab = rep("", ncol(final1[,-c(remove_vec)])),quote = F)
  cat("\n")
  ###eff mat
  eff_mat1=final[,1:(n1+n2)]
  eff_mat2=final[,(1+n1+n2):(n1+n2+n3)]
  gap=abs(ncol(eff_mat2)-ncol(eff_mat1))
  eff_mat1=cbind(eff_mat1,matrix(0,nrow=nrow(eff_mat1),ncol=gap))
  eff_mat12=rbind(eff_mat1,eff_mat2)
  ###Parameters need to add
  list1=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),
          paste("Number of  blocks(b1) =",v),
          paste("Number of sub blocks(b2) =",2*v),
          paste("Block size(k1) =",v-1),paste("1st sub block size(k21) =",(n1+n2)),paste("2nd sub block size(k22) =",n3),
          paste("1st sub-sub block size(k31) =",n1),paste("2nd sub-sub block size(k32) =",n2),paste("2nd sub-sub block size(k33) =",n3),
          paste("Number of times each pair occur together in blocks (lambda_1) =",(m*n-2)),
                paste("Number of times pair of 1st associates occur together in sub blocks (lambda_21) =",m*(n-2)),
                paste("Number of times pair of 2nd associates occur together in sub blocks (lambda_22) =",n*(m-2)),
                paste("Number of times pair of 3rd associates occur together in sub blocks (lambda_23) =",m*(n-2))
          )
  print(list1,quote=F)
  IBD_Cannonical(eff_mat12)
  cat("\n")
}
if(D4==T){
#print("Nested Incomplete Block Design",quote=F)
  message("IBD at block level")
  cat("\n")
print(final)
cat("\n")
list2=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),
        paste("Number of  blocks(b1) =",v),
        paste("Block size(k1) =",v-1),
        paste("Number of times each pair occur together in blocks (lambda_1) =",(m*n-2))
)
print(list2,quote=F)
IBD_Cannonical(final)
cat("\n")
}
if(D5==T){
  message("IBD at sub-block level")
  cat("\n")
  ####eff mat
  eff_mat1=final[,1:n1]
  eff_mat2=final[,(1+n1):(n1+n2)]
  eff_mat3=final[,(1+n1+n2):(n1+n2+n3)]
  gap1=ncol(eff_mat3)-ncol(eff_mat1)
  gap2=ncol(eff_mat3)-ncol(eff_mat2)
  eff_mat1=cbind(eff_mat1,matrix(0,nrow=nrow(eff_mat1),ncol=gap1))
  eff_mat2=cbind(eff_mat2,matrix(0,nrow=nrow(eff_mat2),ncol=gap2))
  final2=rbind(eff_mat1,eff_mat2,eff_mat3)
  print(final2)
  cat("\n")
  list3=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),
          paste("Number of sub blocks(b2) =",2*v),
         paste("1st sub block size(k21) =",(n1+n2)),paste("2nd sub block size(k22) =",n3),
          paste("Number of times pair of 1st associates occur together in sub blocks (lambda_21) =",m*(n-2)),
          paste("Number of times pair of 2nd associates occur together in sub blocks (lambda_22) =",n*(m-2)),
         paste("Number of times pair of 3rd associates occur together in sub blocks (lambda_23) =",m*(n-2)),
          paste("Number of times pair of 1st associates occur together in sub-sub blocks (lambda_31) =",m*(n-2)),
          paste("Number of times pair of 2nd associates occur together in sub-sub blocks (lambda_32) =",n*(m-2)),
          paste("Number of times pair of 3rd associates occur together in sub-sub blocks (lambda_33) =",(m-2)*(n-2))
  )
  print(list3,quote=F)
  IBD_Cannonical(final2)
  cat("\n")
}
if(D6==T){
  message("IBD at sub-sub block level")
  cat("\n")
c3=final[,((1+n1+n2):ncol(final))]
 c1= cbind(final[,(1+n1):(n1+n2)],matrix(0,nrow = nrow(final),ncol=ncol(c3)-n2))
 c2=cbind(final[,(1:n1)],matrix(0,nrow = nrow(final),ncol=ncol(c3)-n1))
  final3=rbind(c3,c1,c2)
  print(final3)
  cat("\n")
  list4=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),
          paste("Number of sub-sub blocks(b3) =",3*v),
          paste("1st sub-sub block size(k31) =",n1),paste("2nd sub-sub block size(k32) =",n2),paste("2nd sub-sub block size(k33) =",n3),
          paste("Number of times each pair occur together in blocks (lambda_1) =",(m*n-2)),
          paste("Number of times pair of 1st associates occur together in sub-sub blocks (lambda_31) =",m*(n-2)),
          paste("Number of times pair of 2nd associates occur together in sub-sub blocks (lambda_32) =",n*(m-2)),
          paste("Number of times pair of 3rd associates occur together in sub-sub blocks (lambda_33) =",(m-2)*(n-2))
  )
  print(list4,quote=F)
  IBD_Cannonical(final3)
}
}

}else{
  message("Please enter m and n such that v = m*n and v>=12")
}
}
#Series3(4,4)

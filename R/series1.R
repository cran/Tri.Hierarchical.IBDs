
#' Tri-Hierarchical IBDs using Triangular Association Scheme
#'@name Series1
#' @param v Number of treatments, v = nC2 where n>=5
#' @param D1 Bi-Hierarchical IBD by ignoring blocks
#' @param D2 Bi-Hierarchical IBD by ignoring sub-blocks
#' @param D3 Bi-Hierarchical IBD by ignoring sub-sub blocks
#' @param D4 IBD at block level
#' @param D5 IBD at sub block level
#' @param D6 IBD at sub-sub block level
#' @param Randomization Randomization of layout of the designs if needed enter TRUE; by default it is FALSE.
#' @description This function generates Tri-Hierarchical IBDs based on Triangular association scheme. Here, v=nC2, n >=5. We find balanced incomplete block designs (BIBD) at block level and triangular PBIB designs at sub-block level as well as sub-sub block level. Information matrix pertaining to the estimation of treatments effects, canonical efficiency factor in comparison to an orthogonal design and six component designs are provided.
#' @return It gives Tri-HIB design and six component designs with canonical efficiency factor in comparison to an orthogonal design.
#' @references
#'Preece, D.A. (1967) <https://doi.org/10.1093/biomet/54.3-4.479>."Nested balanced incomplete block designs".
#' @export
#'@note Numbers in the outer most parentheses represents as block elements, second level parentheses as sub block elements and inner most parentheses as sub-sub block elements.
#'
#' @examples
#' library(Tri.Hierarchical.IBDs)
#' Series1(15,D1=TRUE,D2=TRUE,D3=TRUE,D4=TRUE,D5=FALSE,D6=TRUE,Randomization=FALSE)
Series1=function(v,D1=FALSE,D2=FALSE,D3=FALSE,D4=FALSE,D5=FALSE,D6=FALSE,Randomization=FALSE){
  n=0.5+(sqrt(2*v+0.25))
  if(as.integer(n)-n==0){
  mat=matrix(1:(n)^2,nrow=n,ncol=n)
  mat[lower.tri(mat)]<-0
  mat[diag(mat)]<-0
  ###############
  vec=c(1:v)
  k=1
  for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
      if(mat[i,j]!=0){
        mat[i,j]=k
        k=k+1
      }
    }
  }
  ###############

  for(i in 1:ncol(mat)){
   mat[,i]<-mat[i,]
  }
  ####################

##################
  newmat=NULL
  for(i in 1:v){
  posmat=which(t(mat)==i,arr.ind = T)[1,]
  AS2=c(setdiff(mat[,posmat[1]],i))
  AS2=AS2[AS2>0]
  AS1=c(setdiff(mat[posmat[2],],i))
  AS1=AS1[AS1>0]
  AS3=setdiff(1:v,(c(AS1,AS2,i)))
  AS3
  x=c(AS1,AS2,AS3)
  newmat=rbind(newmat,x)
  }
  row.names(newmat)<-NULL
  ################
  #write.table(newmat,"clipboard-128",sep="\t")
  ###################
  N_N_prime=function(matrix){
    ysd=matrix
  delprime=matrix(0,nrow=length(ysd),ncol=max(ysd))
  for(j in 1:length(ysd)){
    delprime[j,(t(ysd))[j]]<-1
  }
 #########################
  D1_mat_prime=matrix(0,nrow=length(ysd),ncol=nrow(ysd))
  k=1
  for(j in 1:nrow(ysd)){
    D1_mat_prime[(k):(k-1+ncol(ysd)),j]=1
    k=k+ncol(ysd)
  }
  ##########
  N=t(delprime)%*%D1_mat_prime
  NN_p=N%*%t(N)
  return(NN_p)
  }
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
 final=newmat
  n1=length(AS1)
  n2=length(AS2)
  n3=length(AS3)
  n11=n1+n2
  n22=n3
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
  r=v-1
  #N_N_prime(final)
  cat("\n")

  list1=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),
          paste("Number of  blocks(b1) =",v),
          paste("Number of sub blocks(b2) =",2*v), paste("Number of sub-sub blocks(b3) =",3*v),
          paste("Block size(k1) =",v-1),paste("1st sub block size(k21) =",(n11)),paste("2nd sub block size(k22) =",n22),
          paste("1st sub-sub block size(k31) =",n11/2),paste("2nd sub-sub block size(k32) =",n22),
          paste("Number of times each pair occur together in blocks (lambda_1) =",(v-2)),
          paste("Number of times pair of 1st associates occur together in sub blocks (lambda_21) =",v-(n11)),
          paste("Number of times pair of 2nd associates occur together in sub blocks (lambda_22) =",((((n-4)*(n-5))/2)+4)),
          paste("Number of times pair of 1st associates occur together in sub-sub blocks (lambda_31) =",v-(n11+1)),
          paste("Number of times pair of 2nd associates occur together in sub-sub blocks (lambda_32) =",((n-4)*(n-5))/2)
  )
  print(list1,quote=F)
  #####################sub-sub block efficiency
  final_for_eff=rbind(final[,1:n1],final[,(n1+1):(n1+n2)])
  final_for_eff=cbind(final_for_eff,matrix(0,nrow=nrow(final_for_eff),ncol=(n3-n2)))
  final_for_eff=rbind(final_for_eff,final[,(n1+n2+1):(n1+n2+n3)])
  ################
  IBD_Cannonical(final_for_eff)
  cat("\n")
  if(D1==T){
    message("Bi-Hierarchical IBD by ignoring blocks")
    cat("\n")
    #message("Bi-Hierarchical IBD",quote=F)
    final4=final1[,-c(1,ncol(final1))]
    final41=final4[,(1:(n1+n2+4+2))]
    final42=final4[,((1+n1+n2+4+2):ncol(final4))]
    gap=abs(ncol(final41)-ncol(final42))
    gap_mat=matrix(0,nrow=nrow(final),ncol=gap)
    if(ncol(final41)>ncol(final42)){
      final42=cbind(final42,gap_mat)
    }else{
      final41=cbind(final41,gap_mat)
    }
    final4=rbind(final41,final42)


    #,matrix(0,nrow=nrow(final1),ncol=(ncol(final4[,(1:(n1+n2+4+2))])-ncol(final4[,((1+n1+n2+4+2):ncol(final4))])))))

    final4[final4==0]=NA
    blk5=c()
    for(i in 1:nrow(final4)){
      blk5=c(blk5,paste0("Block-",i,collapse = ""))
    }
    row.names(final4)<-blk5
    prmatrix(final4, rowlab = , collab = rep("", ncol(final4)),quote = F,na.print = "")
    cat("\n")
    list5=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),
            paste("Number of sub blocks(b2) =",2*v), paste("Number of sub-sub blocks(b3) =",3*v),
            paste("1st sub block size(k21) =",(n11)),paste("2nd sub block size(k22) =",n22),
            paste("1st sub-sub block size(k31) =",n11/2),paste("2nd sub-sub block size(k32) =",n22),
            paste("Number of times pair of 1st associates occur together in sub blocks (lambda_21) =",v-(n11)),
            paste("Number of times pair of 2nd associates occur together in sub blocks (lambda_22) =",((((n-4)*(n-5))/2)+4)),
            paste("Number of times pair of 1st associates occur together in sub-sub blocks (lambda_31) =",v-(n11+1)),
            paste("Number of times pair of 2nd associates occur together in sub-sub blocks (lambda_32) =",((n-4)*(n-5))/2)
    )

    print(list5,quote=F)
    #############eff mat
    final_eff1=final[,1:n1]
    final_eff2=final[,(1+n1):(n1+n2)]
    final_eff3=final[,(n11+1):(n11+n22)]
    gap=ncol(final_eff3)-ncol(final_eff1)
    final_eff1=cbind(final_eff1,matrix(0,nrow=nrow(final_eff1),ncol=gap))
    final_eff2=cbind(final_eff2,matrix(0,nrow=nrow(final_eff2),ncol=gap))
    final2=rbind(final_eff1,final_eff2,final_eff3)

    IBD_Cannonical(final2)
    cat("\n")

  }
  if(D2==T){
    message("Bi-Hierarchical IBD by ignoring sub blocks")
    cat("\n")
    remove_vec=c(2,(3+n1+2+n2+1+1),(3+n1+2+n2+1+2),(ncol(final1)-1))
    final5=final1[,-c(remove_vec)]
    prmatrix(final5, rowlab = , collab = rep("", ncol(final5)),quote = F,na.message = "")
    cat("\n")

    list5=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),
            paste("Number of  blocks(b1) =",v), paste("Number of sub-sub blocks(b3) =",3*v),
            paste("Block size(k1) =",v-1),paste("1st sub-sub block size(k31) =",n11/2),paste("2nd sub-sub block size(k32) =",n22),
            paste("Number of times each pair occur together in blocks (lambda_1) =",(v-2)),
            paste("Number of times pair of 1st associates occur together in sub-sub blocks (lambda_31) =",v-(n11+1)),
            paste("Number of times pair of 2nd associates occur together in sub-sub blocks (lambda_32) =",((n-4)*(n-5))/2)
    )
    print(list5,quote=F)
    ####eff matrix
    eff_mat1=rbind(final[,1:n1],final[,((n1+1):(n1+n2))])
    eff_mat2=final[,(n1+n2+1):ncol(final)]
    gap=ncol(eff_mat2)-ncol(eff_mat1)
    eff_mat1=cbind(eff_mat1,matrix(0,nrow=nrow(eff_mat1),ncol=gap))
    final5=rbind(eff_mat1,eff_mat2)
    IBD_Cannonical(final5)
    cat("\n")
  }
  if(D3==T){
    message("Bi-Hierarchical IBD by ignoring sub-sub blocks")
    cat("\n")
    remove_vec2=c(3,(3+n1+1),(3+n1+2),(3+n1+2+n2+1),(3+n1+2+n2+1+3),(ncol(final1)-2))
    final6=final1[,-c(remove_vec2)]
    prmatrix(final6, rowlab = , collab = rep("", ncol(final6)),quote = F,na.message = "")
    cat("\n")
    list6=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),
            paste("Number of  blocks(b1) =",v),
            paste("Number of sub blocks(b2) =",2*v),
            paste("Block size(k1) =",v-1),paste("1st sub block size(k21) =",(n11)),paste("2nd sub block size(k22) =",n22),
            paste("Number of times each pair occur together in blocks (lambda_1) =",(v-2)),
            paste("Number of times pair of 1st associates occur together in sub blocks (lambda_21) =",v-(n11)),
            paste("Number of times pair of 2nd associates occur together in sub blocks (lambda_22) =",((((n-4)*(n-5))/2)+4))
    )

    print(list6,quote=F)
    ####eff matrix
    eff_mat1=final[,1:(n1+n2)]
    eff_mat2=final[,((n1+n2+1):(n1+n2+n3))]
    gap=(ncol(eff_mat2)-ncol(eff_mat1))
    if(gap<0){
      a=abs(gap)
      eff_mat2=cbind(eff_mat2,matrix(0,nrow=nrow(eff_mat2),ncol=a))
    }else{
      a=abs(gap)
      eff_mat1=cbind(eff_mat1,matrix(0,nrow=nrow(eff_mat1),ncol=a))
    }

    final6=rbind(eff_mat1,eff_mat2)
    IBD_Cannonical(final6)
    cat("\n")
  }
  if(D4==T){
    message("IBD at block level")
    cat("\n")
    blk1=c()
    for(i in 1:nrow(final)){
      blk1=c(blk1,paste0("Block-",i,collapse = ""))
    }
    row.names(final)<-blk1
    print(final)
    cat("\n")

    list2=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),
            paste("Number of  blocks(b1) =",v),
            paste("Number of sub blocks(b2) =",2*v), paste("Number of sub-sub blocks(b3) =",3*v),
            paste("Block size(k1) =",v-1),paste("Number of times each pair occur together in blocks (lambda_1) =",(v-2))
    )

    print(list2,quote=F)
    IBD_Cannonical(final)
    cat("\n")
  }
  if(D5==T){
    message("IBD at sub-block level")
    cat("\n")
    final21=final[,1:(n1+n2)]
    final22=final[,((1+n1+n2):ncol(final))]
    gap=abs(ncol(final21)-ncol(final22))
    gap_mat=matrix(0,nrow=nrow(final),ncol=gap)
    if(ncol(final21)>ncol(final22)){
      final22=cbind(final22,gap_mat)
    }else{
      final21=cbind(final21,gap_mat)
    }
    final2=rbind(final21,final22)
    blk3=c()
    for(i in 1:nrow(final2)){
      blk3=c(blk3,paste0("Block-",i,collapse = ""))
    }
    row.names(final2)<-blk3
    print(final2)
    #N_N_prime(final2)
    cat("\n")
    list3=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),
            paste("Number of sub blocks(b2) =",2*v),
            paste("1st sub block size(k21) =",(n11)),paste("2nd sub block size(k22) =",n22),
            paste("1st sub-sub block size(k31) =",n11/2),paste("2nd sub-sub block size(k32) =",n22),
            paste("Number of times each pair occur together in blocks (lambda_1) =",(v-2)),
            paste("Number of times pair of 1st associates occur together in sub blocks (lambda_21) =",v-(n11)),
            paste("Number of times pair of 2nd associates occur together in sub blocks (lambda_22) =",((((n-4)*(n-5))/2)+4))

    )
    print(list3,quote=F)
    IBD_Cannonical(final2)
  cat("\n")
  }

  if(D6==T){

    message("IBD at sub-sub block level")
    cat("\n")
    c3=final[,((1+n1+n2):ncol(final))]
    c1= cbind(final[,1:n1],matrix(0,nrow = nrow(final),ncol=ncol(c3)-n1))
    c2=cbind(final[,(1+n1):(n1+n2)],matrix(0,nrow = nrow(final),ncol=ncol(c3)-n2))
    final3=rbind(c3,c1,c2)
    blk4=c()
    for(i in 1:nrow(final3)){
      blk4=c(blk4,paste0("Block-",i,collapse = ""))
    }
    row.names(final3)<-blk4
    print(final3)
    cat("\n")
    list4=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),
            paste("Number of  blocks(b1) =",v),
            paste("Number of sub-sub blocks(b3) =",3*v),
            paste("1st sub-sub block size(k31) =",n11/2),paste("2nd sub-sub block size(k32) =",n22),
            paste("Number of times each pair occur together in blocks (lambda_1) =",(v-2)),
            paste("Number of times pair of 1st associates occur together in sub-sub blocks (lambda_31) =",v-(n11+1)),
            paste("Number of times pair of 2nd associates occur together in sub-sub blocks (lambda_32) =",((n-4)*(n-5))/2)
    )

    print(list4,quote=F)
    IBD_Cannonical(final3)
  }
  }else{
  message("Please enter v such that v = nC2 and n>=5")
  }
}



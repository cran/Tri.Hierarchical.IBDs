
#' Tri-Hierarchical IBDs using Initial Block Solution
#'
#' @param v Number of treatments, (11 <= v < 200) a prime number
#' @param D1 Bi-Hierarchical IBD by ignoring blocks
#' @param D2 Bi-Hierarchical IBD by ignoring sub-blocks
#' @param D3 Bi-Hierarchical IBD by ignoring sub-sub blocks
#' @param D4 IBD at block level
#' @param D5 IBD at sub block level
#' @param D6 IBD at sub-sub block level
#' @param Randomization Randomization of layout of the designs if needed enter TRUE; by default it is FALSE.
#'@description
#'This function gives Tri-Hierarchical IBDs using initial sequences. Here,  v = 4t+1 or 4t+3, where t is an integer and v should be a prime number, using primitive element of Galois field designs are generated. We find balanced incomplete block designs (BIBD) at block level and PBIB designs at sub-block level as well as sub-sub block level with circular association scheme. Information matrix pertaining to the estimation of treatments effects, canonical efficiency factor in comparison to an orthogonal design and six component designs are provided.
#' @references
#'Preece, D.A. (1967) <https://doi.org/10.1093/biomet/54.3-4.479>."Nested balanced incomplete block designs".
#' @return It gives Tri-HIB design and six component designs with canonical efficiency factor in comparison to an orthogonal design.
#' @export
#'@note Numbers in the outer most parentheses represents as block elements, second level parentheses as sub block elements and inner most parentheses as sub-sub block elements.
#'
#' @examples
#' library(Tri.Hierarchical.IBDs)
#' Series4(13,D1=FALSE,D2=FALSE,D3=TRUE,D4=TRUE,D5=FALSE,D6=TRUE,Randomization=TRUE)
Series4=function(v,D1=FALSE,D2=FALSE,D3=FALSE, D4=FALSE,D5=FALSE,D6=FALSE,Randomization=FALSE){

############
  is.prime=function(v){
    if(v==2){
      return(TRUE)
    }
    if(v>2){
      a=c(2:as.integer(v/2))
      b=rep(v,time=length(a))
      c=c(b%%a)
      d=length(c[c==0])

      if(d==0){
        return(TRUE)
      }else{
        return(FALSE)
      }
    }
  }
  if(is.prime(v)==T && v>10){
###############################
t1=((v-3)/4)
t2=((v-1)/4)

#############################
if(t1%%1==0){
 t=t1
}else{
 t=t2
}
#######################
if(t==t1){
  #4t+3
  boolean=T
}else{
  #4t+1
  boolean=F
}
###############################
if(v==3){
  pr=2
}
if(v==5){
  pr=2
}
if(v==7){
  pr=3
}
if(v==11){
  pr=2
}
if(v==13){
  pr=2
}
if(v==17){
  pr=3
}
if(v==19){
  pr=2
}
if(v==23){
  pr=5
}
if(v==29){
  pr=3
}
if(v==31){
  pr=3
}
if(v==37){
  pr=3
}
if(v==41){
  pr=6
}
if(v==43){
  pr=3
}
if(v==47){
  pr=5
}
if(v==53){
  pr=2
}
if(v==59){
  pr=2
}
if(v==61){
  pr=2
}

if(v==67){
  pr=2
}
if(v==71){
  pr=7
}
if(v==73){
  pr=5
}
if(v==79){
  pr=3
}
if(v==83){
  pr=2
}
if(v==89){
  pr=3
}
if(v==97){
  pr=5
}
if(v==101){
  pr=2
}
if(v==103){
  pr=5
}
if(v==107){
  pr=2
}
if(v==109){
  pr=6
}
if(v==113){
  pr=3
}
if(v==127){
  pr=3
}
if(v==131){
  pr=2
}
if(v==137){
  pr=3
}
if(v==139){
  pr=2
}
if(v==149){
  pr=2
}
if(v==151){
  pr=6
}
if(v==157){
  pr=2
}
if(v==163){
  pr=2
}
if(v==167){
  pr=2
}
if(v==173){
  pr=5
}
if(v==179){
  pr=2
}
if(v==181){
  pr=2
}
if(v==191){
  pr=19
}
if(v==193){
  pr=5
}
if(v==197){
  pr=2
}
if(v==199){
  pr=7
}
###########################
if(boolean==F){
  #4t+1
even=seq(0,(v-3),by=2)
odd=seq(1,(v-2),by=2)
}else{
  #4t+3
  even=seq(0,(v-1),by=2)
  even=even[-c(length(even))]
  odd=seq(1,(v-2),by=2)
}
####################################################

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

####################%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vec=NULL
for(i in odd){
  vec=c(vec,pr^i)
}
vec=vec%%v
first_part=NULL
for(i in 0:(v-1)){
  first_part=rbind(first_part,c(vec+i))
}
first_part=first_part%%v
first_part[first_part==0]<-v
###################################
second_part=NULL
vec1=NULL
for(i in even){
  vec1=c(vec1,pr^i)
}
vec1=vec1%%v
second_part=NULL
for(i in 0:(v-1)){
  second_part=rbind(second_part,c(vec1+i))
}
second_part=second_part%%v
second_part[second_part==0]<-v
#############################
design_final=cbind(second_part,first_part)
DD1=design_final
#####################Randomization part
rand=c(sample(1:nrow(design_final)))
if(Randomization==T){
  design_final=design_final[c(rand),]
}
###########################$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
D11=DD1
#write.table(D11,"clipboard",sep="\t")
#####################################we need earlier so we need to add here
last_design=rbind(D11[,((t+(ncol(second_part)-t))+1):ncol(D11)],cbind(D11[,1:t],matrix(0,nrow=nrow(D11),ncol=ncol(D11[,((t+(ncol(second_part)-t))+1):ncol(D11)])-
                                                                   ncol(D11[,1:t]))),
                  cbind(D11[,(1+t):ncol(second_part)],matrix(0,nrow=nrow(D11),ncol=ncol(D11[,((t+(ncol(second_part)-t))+1):ncol(D11)])-
                                                               ncol(D11[,(1+t):ncol(second_part)]))))
###########incidence matrix
incidence=matrix(0,nrow=v,ncol=nrow(last_design))
for(i in 1:v){
  pos=which(last_design==i,arr.ind = T)[,1]
  incidence[i,c(pos)]<-1
}

##############################
N=incidence
NNp=N%*%t(N)
#############################
lambdas=NNp[1,2:(((v-1)/2)+1)]
###################
array_lambda=NULL
for(i in 1:length(lambdas)){
  array_lambda=rbind(array_lambda,paste0("lambda_3",i," = ",lambdas[i]))
}
array_lambda=noquote(array_lambda)

#####we need earlier so we need to add here


#########################$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
n1=length(1:t)
n2=length((t+1):ncol(second_part))
n3=length(((t+(ncol(second_part)-t))+1):ncol(design_final))
n11=n1+n2
n22=n3

final=design_final
design_final=cbind("{","[","(",design_final[,1:t],")","(", design_final[,(1+t):ncol(second_part)],")","]",
                   "[","(",design_final[,((t+(ncol(second_part)-t))+1):ncol(design_final)],")","]","}")
message("Tri-Hierarchical IBD")
cat("\n")
prmatrix(design_final, rowlab = , collab = rep("", ncol(design_final)),quote = F)
cat("\n")
if(boolean==T){
  #4t+3
list1=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),
        paste("Number of  blocks(b1) =",v),
        paste("Number of sub blocks(b2) =",2*v), paste("Number of sub-sub blocks(b3) =",3*v),
        paste("Block size(k1) =",v-1),paste("Sub block size(k2) =",(2*t)+1),
        paste("1st sub-sub block size(k31) =",t),paste("2nd sub-sub block size(k32) =",t+1),
        paste("3rd sub-sub block size(k33) =",(2*t)+1),
        paste("Number of times each pair occur together in blocks (lambda_1) =",(v-2)),
        paste("Number of times each pair occur together in  sub blocks (lambda_2) =",(2*t)),
        array_lambda
)

}else{
  #4t+1
  list1=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),
          paste("Number of  blocks(b1) =",v),
          paste("Number of sub blocks(b2) =",2*v), paste("Number of sub-sub blocks(b3) =",3*v),
          paste("Block size(k1) =",v-1),paste("Sub block size(k2) =",2*t),
          paste("1st sub-sub block size(k31) =",t),paste("2nd sub-sub block size(k32) =",t),
          paste("3rd sub-sub block size(k33) =",2*t),
          paste("Number of times each pair occur together in blocks (lambda_1) =",(v-2)),
          paste("Number of times each pair occur together in  sub blocks (lambda_2) =",(2*t-1)),
          array_lambda
  )
}
print(list1,quote=F)
s1=second_part[,1:t]
s2=second_part[,(t+1):ncol(second_part)]
s3=final[,(1+ncol(second_part)):ncol(final)]
gap1=ncol(s3)-ncol(s1)
gap2=ncol(s3)-ncol(s2)
mat_eff1=cbind(s1,matrix(0,nrow=nrow(s1),ncol=gap1))
mat_eff2=cbind(s2,matrix(0,nrow=nrow(s2),ncol=gap2))
IBD_Cannonical(rbind(mat_eff1,mat_eff2,s3))
cat("\n")

##################################################
# print("Design 2",quote=F)
# D21=D2[,1:(2*t)]
# D22=D2[,((2*t)+1):ncol(design_final)]
# D21=cbind("[","(",D21[,1:ncol(D21)],")","]")
# gap=7+2*t
if(D1==T){
  message("Bi-Hierarchical IBD by ignoring blocks")
  cat("\n")
  remove_vec=c(1,ncol(design_final))
  prmatrix(design_final[,-c(remove_vec)], rowlab = , collab = rep("", ncol(design_final[,-c(remove_vec)])),quote = F)
  ##########
  cat("\n")
  if(boolean==T){
    #4t+3
    list1=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),
            paste("Number of sub blocks(b2) =",2*v), paste("Number of sub-sub blocks(b3) =",3*v),
            paste("Sub block size(k2) =",(2*t)+1),
            paste("1st sub-sub block size(k31) =",t),paste("2nd sub-sub block size(k32) =",t+1),
            paste("Number of times each pair occur together in  sub blocks (lambda_2) =",(2*t)),
            array_lambda
    )

  }else{
    #4t+1
    list1=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),
            paste("Number of sub blocks(b2) =",2*v), paste("Number of sub-sub blocks(b3) =",3*v),
            paste("Sub block size(k2) =",2*t),
            paste("1st sub-sub block size(k31) =",t),paste("2nd sub-sub block size(k32) =",t),
            paste("Number of times each pair occur together in  sub blocks (lambda_2) =",(2*t-1)),
            array_lambda
    )
  }
  print(list1,quote=F)
  ####eff matrix
  eff_mat1=final[,1:n1]
  eff_mat2=final[,(1+n1):(n1+n2)]
  eff_mat3=final[,(n1+n2+1):(n1+n2+n3)]
  gap1=ncol(eff_mat3)-ncol(eff_mat1)
  gap2=ncol(eff_mat3)-ncol(eff_mat2)
  eff_mat11=cbind(eff_mat1,matrix(0,nrow=nrow(eff_mat1),ncol=gap1))
  eff_mat22=cbind(eff_mat2,matrix(0,nrow=nrow(eff_mat2),ncol=gap2))
  eff_mat=rbind(eff_mat11,eff_mat22,eff_mat3)
  IBD_Cannonical(eff_mat)
  cat("\n")
}
if(D2==T){
  message("Bi-Hierarchical IBD by ignoring sub blocks")
  cat("\n")
  remove_vec=c(2,(3+n1+2+n2+1+1),(3+n1+2+n2+1+2),ncol(design_final)-1)
  prmatrix(design_final[,-c(remove_vec)], rowlab = , collab = rep("", ncol(design_final[,-c(remove_vec)])),quote = F)
  ##########
  cat("\n")
  if(boolean==T){
    #4t+3
    list1=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),
            paste("Number of  blocks(b1) =",v),
             paste("Number of sub-sub blocks(b3) =",3*v),
            paste("Block size(k1) =",v-1),
            paste("1st sub-sub block size(k31) =",t),paste("2nd sub-sub block size(k32) =",t+1),
            paste("Number of times each pair occur together in blocks (lambda_1) =",(v-2)),
            array_lambda
    )

  }else{
    #4t+1
    list1=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),
            paste("Number of  blocks(b1) =",v),
            paste("Number of sub-sub blocks(b3) =",3*v),
            paste("Block size(k1) =",v-1),
            paste("1st sub-sub block size(k31) =",t),paste("2nd sub-sub block size(k32) =",t),
            paste("Number of times each pair occur together in blocks (lambda_1) =",(v-2)),
            array_lambda
    )
  }
  print(list1,quote=F)
  ####eff matrix
  eff_mat1=final[,1:n1]
  eff_mat2=final[,(1+n1):(n1+n2)]
  eff_mat3=final[,(n1+n2+1):(n1+n2+n3)]
  gap1=ncol(eff_mat3)-ncol(eff_mat1)
  gap2=ncol(eff_mat3)-ncol(eff_mat2)

  eff_mat11=cbind(eff_mat1,matrix(0,nrow=nrow(eff_mat1),ncol=gap1))
  eff_mat22=cbind(eff_mat2,matrix(0,nrow=nrow(eff_mat2),ncol=gap2))
  eff_mat=rbind(eff_mat11,eff_mat22,eff_mat3)
  IBD_Cannonical(eff_mat)
  cat("\n")
}
if(D3==T){
  message("Bi-Hierarchical IBD by ignoring sub-sub blocks")

  cat("\n")
  remove_vec=c(3,(3+n1+1),(3+n1+2),(3+n1+2+n2+1),(3+n1+2+n2+1+2+1),ncol(design_final)-2)
  prmatrix(design_final[,-c(remove_vec)], rowlab = , collab = rep("", ncol(design_final[,-c(remove_vec)])),quote = F)
  cat("\n")
  ##########
  if(boolean==T){
    #4t+3
    list1=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),
            paste("Number of  blocks(b1) =",v),
            paste("Number of sub blocks(b2) =",2*v),
            paste("Block size(k1) =",v-1),paste("Sub block size(k2) =",(2*t)+1),
            paste("Number of times each pair occur together in blocks (lambda_1) =",(v-2)),
            paste("Number of times each pair occur together in  sub blocks (lambda_2) =",(2*t))

    )

  }else{
    #4t+1
    list1=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),
            paste("Number of  blocks(b1) =",v),
            paste("Number of sub blocks(b2) =",2*v),
            paste("Block size(k1) =",v-1),paste("Sub block size(k2) =",2*t),
            paste("Number of times each pair occur together in blocks (lambda_1) =",(v-2)),
            paste("Number of times each pair occur together in  sub blocks (lambda_2) =",(2*t-1))

    )
  }
  print(list1,quote=F)
  IBD_Cannonical(rbind(final[,1:n11],final[,(n11+1):(n11+n22)]))
  cat("\n")
}
if(D4==T){
  message("IBD at block level")
cat("\n")
##################################################################
#DD1=cbind("{",DD1,"}")
#prmatrix(DD1, rowlab = , collab = rep("", ncol(DD1)),quote = F,na.print = "")
print(D11)
cat("\n")
################################
if(boolean==T){
  #4t+3
  list1=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),
          paste("Number of  blocks(b1) =",v),
          paste("Block size(k1) =",v-1),
          paste("Number of times each pair occur together in blocks (lambda_1) =",(v-2))
  )

}else{
  #4t+1
  list1=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),
          paste("Number of  blocks(b1) =",v),
          paste("Block size(k1) =",v-1),
          paste("Number of times each pair occur together in blocks (lambda_1) =",(v-2))

  )
}

##################################
print(list1,quote=F)
IBD_Cannonical(D11)
cat("\n")
}
if(D5==T){
  message("IBD at sub-block level")
  cat("\n")
  design1=final[,1:(n1+n2)]
  design2=final[,(n1+n2+1):(n1+n2+n3)]
  gap=ncol(design2)-ncol(design1)
  design1=cbind(design1,matrix(0,nrow=nrow(design1),ncol=gap))
print(rbind(design1,design2))
cat("\n")
#################################
if(boolean==T){
  #4t+3
  list1=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),
          paste("Number of sub blocks(b2) =",2*v), paste("Number of sub-sub blocks(b3) =",3*v),
         paste("Sub block size(k2) =",(2*t)+1),
          paste("1st sub-sub block size(k31) =",t),paste("2nd sub-sub block size(k32) =",t+1),
          paste("3rd sub-sub block size(k33) =",(2*t)+1),

          paste("Number of times each pair occur together in  sub blocks (lambda_2) =",(2*t))

  )

}else{
  #4t+1
  list1=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),

          paste("Number of sub blocks(b2) =",2*v), paste("Number of sub-sub blocks(b3) =",3*v),
          paste("Sub block size(k2) =",2*t),
          paste("1st sub-sub block size(k31) =",t),paste("2nd sub-sub block size(k32) =",t),
          paste("3rd sub-sub block size(k33) =",2*t),

          paste("Number of times each pair occur together in  sub blocks (lambda_2) =",(2*t-1))

  )
}

######################
print(list1,quote=F)
IBD_Cannonical(rbind(design1,design2))
cat("\n")
}
#############################
if(D6==T){

  message("IBD at sub-sub block level")
  cat("\n")
#first_part_1=cbind("[",first_part,"]")
#second_part_1=cbind("[",second_part[,1:t],second_part[,(t+1):ncol(second_part)],"]")
  deg1=final[,1:n1]
  deg2=final[,(1+n1):(n1+n2)]
  deg3=final[,(n1+n2+1):(n1+n2+n3)]
gap1=ncol(deg3)-ncol(deg1)
gap2=ncol(deg3)-ncol(deg2)
deg1=cbind(deg1,matrix(0,nrow=nrow(deg1),ncol=gap1))
deg2=cbind(deg2,matrix(0,nrow=nrow(deg2),ncol=gap2))
deg=rbind(deg1,deg2,deg3)
print(deg)
cat("\n")
#prmatrix(final_D4, rowlab = , collab = rep("", ncol(final_D4)),quote = F,na.print = "")

#################################
if(boolean==T){
  #4t+3
  list1=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),

          paste("Number of sub blocks(b2) =",2*v),
         paste("Sub block size(k2) =",(2*t)+1),
          paste("Number of times each pair occur together in  sub blocks (lambda_2) =",(2*t))
  )

}else{
  #4t+1
  list1=c(paste("Number of treatments(v) =",v), paste("Number of replications(r) =",(v-1)),

          paste("Number of sub blocks(b2) =",2*v),
          paste("Sub block size(k2) =",2*t),
          paste("Number of times each pair occur together in  sub blocks (lambda_2) =",(2*t-1))

  )

}


##########################
print(list1,quote=F)
IBD_Cannonical(deg)
cat("\n")
}
##################################
  }else{
  message("Please enter a prime number, v (<200) such that v = 4t+1 or 4t+3, where t>=2")
  }
}



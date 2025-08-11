#########################################################
#### gene name to gene id && ensembl name to gene id ####
#########################################################

suppressMessages({
library(data.table);
library(dplyr);
library(tidyr)
})


#############################################
####     common functions                ####
#############################################

# SMR analysis between two molecular traits
rm_mhc1=function(smr,mhcStart=25000000,mhcEnd=36000000){
    idx=which(smr$Expo_Chr==6 & ((smr$Expo_bp<=mhcEnd & smr$Expo_bp>=mhcStart) | (smr$topSNP_bp<=mhcEnd & smr$topSNP_bp>=mhcStart)))
    if(length(idx)>0){
        smr=smr[-idx,]
    }
    idx=which(smr$Outco_Chr==6 & smr$Outco_bp<=mhcEnd & smr$Outco_bp>=mhcStart)
    if(length(idx)>0){
        smr=smr[-idx,]
    }
    return(smr)
}

# smr analysis in trans-region
rm_mhc2=function(smr,mhcStart=25000000,mhcEnd=36000000){
    idx=which(smr$ProbeChr==6 & (smr$Probe_bp<=mhcEnd & smr$Probe_bp>=mhcStart))
    if(length(idx)>0){
        smr=smr[-idx,]
    }
    idx=which(smr$topSNP_chr==6 & (smr$topSNP_bp<=mhcEnd & smr$topSNP_bp>=mhcStart))
    if(length(idx)>0){
        smr=smr[-idx,]
    }
    idx=which(smr$trans_chr==6 & ((smr$trans_rightBound<=mhcEnd & smr$trans_rightBound>=mhcStart) | (smr$trans_leftBound<=mhcEnd & smr$trans_leftBound>=mhcStart) | (smr$trans_rightBound<=mhcEnd & smr$trans_leftBound>=mhcStart)))
    if(length(idx)>0){
        smr=smr[-idx,]
    }
    return(smr)
}

mergeframe = function(data_a, data_b, range, nm) {
  buf   <- ReOrder(data_b, colnames(data_b)==nm, data_a$nm)
  i_buf <- which(colnames(data_b)==nm)
  return(merge(data_a, data_b[,c(i_buf, range)], by=nm, all=TRUE))
}

Coefflm = function(y,x,n) {
  res  <- summary(lm(y~x))
  b    <- res$coefficient[2,1]
  se   <- res$coefficient[2,2]
  p    <- res$coefficient[2,4]
  ss   <- res$df[2] + n + 1
  return(list(eff=b, se=se, p=p, size=ss))
}

Residue=function(y,x) {
  b=coefficients(summary(lm(y~., data=x)))
  e=(y-as.matrix(x)%*%b[-1,])[,1]
  return(e)
}

ZScore=function(y,x) {
  sd=sd(y, na.rm=T)
  mean=mean(y,na.rm=T)
  b=coefficients(summary(lm(y~., data=x)))
  e=(y-as.matrix(x)%*%b[-1,])[,1]
  return((e-mean(e,na.rm=T))/sd(e,na.rm=T))
}

ReOrder=function(xMat, i, iVec, NaFlag=TRUE) {
  indx  <- match(iVec,xMat[,i])
  if(NaFlag==TRUE) {
    xBuf  <- xMat[indx[which(!is.na(indx))], ]
  } else {
    xBuf  <- xMat[indx, ]
  }
  return(xBuf)
}

CommonItem=function(list_buf) {
  return(Reduce(intersect, list_buf))
}

OutTest<-function(xVec, d) {
  sd <- sd(xVec, na.rm=T)
  mn <- mean(xVec, na.rm=T)
  xVec[which(abs(xVec-mn)>sd*d)]=NA
  return(xVec)
}

InvNorm=function(x) {
  return(qnorm((rank(x, na.last="keep")-0.5)/sum(!is.na(x))))
}

IstCol=function(xMat, SN, xVec) {
  return(cbind(xMat[,1:SN],xVec,xMat[,-c(1:SN)]))
}

IstRow=function(xMat, SN, xVec) {
  return(rbind(xMat[1:SN,],xVec,xMat[-c(1:SN,)]))
}

# t.test = function(est, se, df) {
#   return( pt(abs(est/se), df, lower.tail=F)*2 )
# }

t.itval = function(est, se, df) {
  d_buf <- qt(0.975, df)*se
  return(list(low=est-d_buf, up=est+d_buf))
}

chi.test = function(chi, df){
  return( pchisq(chi, df, lower.tail=F) )
}

lmpvalue   <-  function(smary)
{
  fstat    <-  smary$fstatistic
  p <- pf(fstat[1],fstat[2],fstat[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

gnrt.norm = function(n, mean=0, sd=1) {
  buf  <- rnorm(n, mean=0, sd=1)
  buf  <- (buf - mean(buf)) /  sd(buf)
  return( (buf + mean) * sd )
}

rmvbinom_sigma=function(size, p, cova){
  library(MASS)
  m  = length(p)
  miu = rep(0, m)
  rawvars = mvrnorm(n=size, mu=miu, Sigma=cova)
  pvars = pnorm(rawvars)
  x  = qbinom(pvars, prob=p, size=2)
}

rmvbinom_old=function(size, p, rho){
  library(bindata)
  #bincorr  <-  (1-r)*diag(2)+r
  
  buf      <-   rmvbin(size*2, p, bincorr=rho)
  sn       <-   c(1:(size*2))
  tmp      <-   2*(sn/2 - as.integer(sn/2))
  buf1     <-   buf[which(tmp==1), ]
  buf2     <-   buf[which(tmp==0), ]
  x        <-   buf1 + buf2  
}

rmvbinom_sigma_old=function(size, p, sigma){
  library(bindata)
  
  buf      <-   rmvbin(size*2, p, sigma=sigma)
  sn       <-   c(1:(size*2))
  tmp      <-   2*(sn/2 - as.integer(sn/2))
  buf1     <-   buf[which(tmp==1), ]
  buf2     <-   buf[which(tmp==0), ]
  x        <-   buf1 + buf2
}

rbvbinom=function(n, p1, p2, r){
  if(r>0){
    fi=sqrt(p1*p2/((1-p1)*(1-p2)));
    gamma=r/(r+fi);
    alpha=p1/(1-gamma);
    beta=p2/(1-gamma);
    k=rbinom(n, 2, gamma);
    k2=which(k==2)
    k01=which(k!=2);
    x=y=k;
    x[k01]=rbinom(length(k01), 2-k[k01], alpha);
    x[k2]=0;
    y[k01]=rbinom(length(k01), 2-k[k01], beta);
    y[k2]=0;
    return(data.frame(x,y));
  }
  if(r<0){
    p2b=p2/(1-p2);
    fi=p2b*sqrt(p1*p2/((1-p1)*(1-p2)));
    gamma=r/(r-fi);
    tau=1-p2/(1-gamma);
    theta=p1/(1-gamma);
    #if(theta>1) theta=1/theta;  # added
    omega=p1*tau/(1-p2);
    if(omega<0) omega=0;
    k=rbinom(n, 2, p2);
    k2=which(k==2);
    k1=which(k==1);
    k0=which(k==0);
    x=y=k;
    x[k2]=rbinom(length(k2), 2, theta);
    y[k2]=0;
    x[k1]=rbinom(length(k1), k[k1], theta)+rbinom(length(k1), 2-k[k1], omega);
    y[k1]=rbinom(length(k1), 2-k[k1], p2b);
    x[k0]=rbinom(length(k0), 2, omega);
    y[k0]=rbinom(length(k0), 2, p2b);
    return(data.frame(x,y));
  }
  if(r==0){
    return(data.frame(x=rbinom(n, 2, p1), y=rbinom(n, 2, p2)));
  }
}

StrandFlip = function(ale, sn) {
  seqbuf = c("A", "T", "A", "C", "G", "C")
  ale  = toupper(ale)
  ibuf = as.numeric()
  ibuf[which(ale=="A")] = 1;  ibuf[which(ale=="T")] = 2;
  ibuf[which(ale=="C")] = 4;  ibuf[which(ale=="G")] = 5;
  ibuf[sn] = ibuf[sn] + 1
  ale  = seqbuf[ibuf]
  return(ale)
}

UpdateAllele_old = function(beta, ale, ref_ale, vecFlag=TRUE) {
  # strand issue
  ale = as.character(ale)
  pb_a1 = as.character(ref_ale[,1]); pb_a2 = as.character(ref_ale[,2])
  sn   = which(ale!=pb_a1&ale!=pb_a2)
  if(length(sn)>0) {
    ale = StrandFlip(ale, sn)
  }
  # reference allele
  sn   = which(ale!=pb_a1)
  if(length(sn)>0) {
    if(vecFlag) {
      beta[sn] = -beta[sn]
    } else beta[sn, ] = -beta[sn,]
  }

  return(beta)
}

UpdateAllele = function(beta, ale, ref_ale, snplist, vecFlag=TRUE) {
  # strand issue
  ale[,1] = toupper(as.character(ale[,1]))
  ale[,2] = toupper(as.character(ale[,2]))
  ref_ale[,1] = toupper(as.character(ref_ale[,1]))
  ref_ale[,2] = toupper(as.character(ref_ale[,2]))
  # multiple alleles
  snbuf1 = which( (ale[,1]==ref_ale[,1] & ale[,2]!=ref_ale[,2]) 
                | (ale[,1]!=ref_ale[,1] & ale[,2]==ref_ale[,2]) 
                | (ale[,1]==ref_ale[,2] & ale[,2]!=ref_ale[,1])
                | (ale[,1]!=ref_ale[,2] & ale[,2]==ref_ale[,1]) )
  snbuf2 = which(ale[,1]!=ref_ale[,1] & ale[,1]!=ref_ale[,2])
  snbuf = unique(c(snbuf1, snbuf2))
  if(length(snbuf)) {
    if(vecFlag) {
      beta = beta[-snbuf];
    } else {
      beta = beta[-snbuf,];
    }
    ale = ale[-snbuf,];
    ref_ale = ref_ale[-snbuf,];
    snplist = snplist[-snbuf]
  }
  # flip the strand
  #snbuf = which(ale[,1]!=ref_ale[,1] & ale[,1]!=ref_ale[,2])
  #if(length(snbuf)) {
  #  ale[,1] = StrandFlip(ale[,1], snbuf)
  #  ale[,2] = StrandFlip(ale[,2], snbuf)
  #}
  # flip the allele
  snbuf = which(ale[,1]!=ref_ale[,1])
  if(length(snbuf)) {
    if(vecFlag) {
      beta[snbuf] = -beta[snbuf] 
    } else {
      beta[snbuf,] = -beta[snbuf,]
    }
  }
  return(list(beta = beta, snplist = snplist))
}

#stepwise VIF function used below
vif_func<-function(in_frame,thresh=10,trace=F){
  library(MASS)
  require(fmsb)
  
  if(class(in_frame) != 'maxtrix') in_frame<-as.matrix(in_frame)
  
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  
  inv_dat = ginv(in_frame)
  mat_q = diag(inv_dat %*% in_frame)
  vif_init = diag(inv_dat)
  vif_init[which(abs(mat_q-1) > 1e-5)] = 1e8 
  vif_max  = max(vif_init)
  
  ndim = dim(in_frame)[1]
  remain = c(1:ndim)
  exclude = as.numeric()
  
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    return(remain)
  } else{
    vif_vals = vif_init
    #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max > thresh){
      max_row = which.max(vif_vals)
      exclude = c(exclude, remain[max_row])
      remain = remain[-max_row]
      vif_vals<-NULL
      in_dat = in_frame[remain, remain]
      vif_vals = diag(ginv(in_dat))
      vif_max  = max(vif_vals)
      if(vif_max<thresh) break
      if(trace==T){ #print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
        flush.console()
      }
    }    
    return(remain)
  }
}

movefile <- function(from, to) {
  todir <- dirname(to)
  if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
  file.rename(from = from,  to = to)
}

UpdateSummaryData = function(sdata1, sdata2) {
  # SNPs in common
  # Columns are SNP, A1, A2, FREQ, B, SE, P, N
  # Consider sdata1 as reference data
  colnames(sdata1) = toupper(colnames(sdata1))
  colnames(sdata2) = toupper(colnames(sdata2))
  snpbuf = Reduce(intersect, list(as.character(sdata1$SNP), 
                                  as.character(sdata2$SNP)))
  sdata1 = sdata1[match(snpbuf, as.character(sdata1$SNP)),]
  sdata2 = sdata2[match(snpbuf, as.character(sdata2$SNP)),]
  ale1 = cbind(as.character(sdata1$A1), as.character(sdata1$A2))
  ale2 = cbind(as.character(sdata2$A1), as.character(sdata2$A2))
  b2 = as.numeric(as.character(sdata2$B))
  # do the alignment
  resbuf = UpdateAllele(b2, ale2, ale1, snpbuf)
  b2 = as.numeric(as.character(resbuf$beta))
  snpbuf = as.character(resbuf$snplist)
  sdata1 = sdata1[match(snpbuf, as.character(sdata1$SNP)),]
  sdata2 = sdata2[match(snpbuf, as.character(sdata2$SNP)),]
  ale1 = cbind(as.character(sdata1$A1), as.character(sdata1$A2))
  se2 = as.numeric(as.character(sdata2$SE))
  #zscore2 = b2/se2
  #pval2 = pchisq(zscore2^2, 1, lower.tail=F)
  pval2 = as.numeric(as.character(sdata2$P))
  b2_old = as.numeric(as.character(sdata2$B))
  snbuf = which(b2*b2_old<0)
  freq = as.numeric(as.character(sdata2$FREQ))
  if(length(snbuf)>0) {
    freq[snbuf] = 1 - freq[snbuf]
  }
  sdata2_new = cbind(snpbuf, ale1, freq, 
                     b2, se2, pval2, as.numeric(as.character(sdata2$N)))
  colnames(sdata2_new) = colnames(sdata2)
  return(list(sdata1=sdata1, sdata2=sdata2_new))
}

n2l = function(seq, seq_type) {
  if (!(seq_type %in% c("prot", "dna", "rna"))) 
    stop("The value of 'what' must be: 'dna', 'rna' or 'prot'.")
  names_list = switch(seq_type, rna = c("a", "c", "g", "u"), 
                       dna = c("a", "c", "g", "t"), 
                      prot = c("a", "b", "c", "d", "e", "f", "g",
                               "h", "i", "j", "k", "l", "m", "n",
                               "o", "p", "q", "r", "s", "t", "v", 
                               "w", "y", "z"))
  elements_list = 1L:length(names_list)
  names(elements_list) = names_list
  degenerate(seq, elements_list)
}

degenerate = function(seq, element_groups) {
  tmp_seq <- seq
  if (!all(unique(tmp_seq) %in% unlist(element_groups))) {
    warning("'seq' contains elements not present in 'element_groups'. Such elements will be replaced by NA.")
    tmp_seq[!(tmp_seq %in% unlist(element_groups))] <- NA
  }
  
  if(is.null(names(element_groups))) {
    warning("'element_groups' is unnamed. Assumed names of groups are their ordinal numbers.")
    names(element_groups) <- 1L:length(element_groups)
  }
  
  if(length(unique(names(element_groups))) != length(names(element_groups))) {
    stop("'element_groups' must have unique names.")
  }
  
  for (i in 1L:length(element_groups)) {
    tmp_seq[tmp_seq %in% element_groups[[i]]] <- names(element_groups)[i]
  }
  
  if(class(seq) == "matrix")
    dim(tmp_seq) <- dim(seq)
  
  tmp_seq
}

get_flip_data<-function(data1,data2){
	# CHECK COLNAMES OF DATA
	if(!all(colnames(data1)==c("SNP","A1","A2","EAF","Beta","SE","P","N"))){
		print("please check the format of data1!")
		print("The columns are SNP A1 A2 EAF Beta SE P N")
		return(NA)
	}
	if(!all(colnames(data2)==c("SNP","A1","A2","EAF","Beta","SE","P","N"))){
		print("please check the format of data2")
		print("The columns are SNP A1 A2 EAF Beta SE P N")
		return(NA)
	}
	# GET COMMON SNPS
	com_idx=match(data1$SNP,data2$SNP,nomatch=0)
	if(length(which(com_idx!=0))==0){
		print("There is no common SNPs between two summary data");
        return(NA)
	}
	data1_com=data1[which(com_idx!=0),]
	data2_com=data2[com_idx,]
	# FLIP
	idx1=which(data1_com$A1==data2_com$A1 & data1_com$A2==data2_com$A2)
	idx2=which(data1_com$A2==data2_com$A1 & data1_com$A1==data2_com$A2)
	data1_com_final=data1_com
	data2_com_final=data2_com
	if(length(idx2)>0){
		data2_com_final$A1[idx2]=data2_com$A2[idx2]
		data2_com_final$A2[idx2]=data2_com$A1[idx2]
		data2_com_final$Beta[idx2]=-data2_com$Beta[idx2]
		if(length(which(!is.na(data2_com_final$EAF)))>0){
			data2_com_final$EAF[idx2]=1-data2_com$EAF[idx2]
		}
	}
	data1_com_final=data1_com_final[c(idx1,idx2),]
	data2_com_final=data2_com_final[c(idx1,idx2),]
	# ALLELE FREQUENCE CHECK
	if(length(which(!is.na(data1_com_final$EAF)))>0 & length(which(!is.na(data2_com_final$EAF)))>0){
		index=which(abs(data1_com_final$EAF-data2_com_final$EAF)>0.2)
	        if(length(index)>0){
			data1_QC=data1_com_final[-index,]
			data2_QC=data2_com_final[-index,]
			QC=list(data1=data1_QC,data2=data2_QC)
		}else{
			QC=list(data1=data1_com_final,data2=data2_com_final)
                } 
	}else{
		QC=list(data1=data1_com_final,data2=data2_com_final)
	}
	return(QC)
}

make_corr_plot<-function(r_est,se_r_est,type="full"){
	lab=0.8;cl = 1;txt_cex=0.8;num=0.8
	library(corrplot)
	ntissue1=nrow(r_est)
	ntissue2=ncol(r_est)
	beta = r_est
	mid_r=(min(beta)+1)/2
	low=mid_r-0.05
	beta = r_est
	mid_r=(min(beta)+1)/2
	low=mid_r-0.05
	high=mid_r+0.05

	min = min(beta) * 0.99
	max = max(abs(beta))*1.00
	barlabels=seq(max, min, length=11)
	trans<-(r_est-mean(c(min,max)))*(-1)
	beta = trans
	min = min(beta) * 0.99
	max = max(abs(beta))*1.00
	trans<-as.matrix(trans)
	corrplot(trans, is.corr = FALSE, method = "color", cl.cex=cl, cl.pos="r",
         cl.lim = c(-max, max), cl.lab = barlabels, #cl.lablim=c(exp(max), exp(-max)),
         tl.srt=45, tl.col="dark blue", mar=c(4,2,6,4), tl.cex=lab,
         cl.ratio=0.07, cl.align.text="l",type=type)
	for( i in 1 : (ntissue1) ) {
    		if(type=="upper"){
        j_index=i:ntissue2
    }
    if(type=="full"){
        j_index = 1 : ntissue2
    }
		for( j in j_index) {
    		if(j<=ntissue2 & r_est[i, j]!=1){
    			if(r_est[i, j]-1.96*se_r_est[i, j]<0){
                output1=output2=""
                }else{
                output1 = sprintf("%0.3f", r_est[i, j])
		        output2 = paste("(", sprintf("%0.2f", se_r_est[i, j]),")", sep="")
		        }
                if(r_est[i,j]<high & r_est[i,j]>low) {
		            text( j, ntissue1-i+1+0.25, output1, cex=num, col="black")
		            text( j, ntissue1-i+1-0.25, output2, cex=txt_cex, col="black")
		        } else {
		            text( j, ntissue1-i+1+0.25, output1, cex=num, col="white")
		            text( j, ntissue1-i+1-0.25, output2, cex=txt_cex, col="white")
		        }
    		}
	    }
	}
}

calcu_std_b_se<-function(z,p,n){
	std_b_hat=z/sqrt(2*p*(1-p)*(n+z^2))
	std_se=1/sqrt(2*p*(1-p)*(n+z^2))
	res<-data.frame(std_b_hat,std_se);
	return(res)
}

calcu_cor_true<-function(b1,se1,b2,se2,rp,jackknife=TRUE){
    idx=which(is.infinite(b1) | is.infinite(b2) | is.infinite(se1) | is.infinite(se2))
    if(length(idx)>0){
        b1=b1[-idx];se1=se1[-idx]
        b2=b2[-idx];se2=se2[-idx]
	rp=rp[-idx]
    }

    var_b1=var(b1,na.rm=T)-mean(se1^2,na.rm=T)
    var_b2=var(b2,na.rm=T)-mean(se2^2,na.rm=T)
    if(var_b1<0){
	var_b1=var(b1,na.rm=T)
}
if(var_b2<0){
	var_b2=var(b2,na.rm=T)
}
    cov_b1_b2=cov(b1,b2,use="complete.obs")-mean(rp,na.rm=T)*sqrt(mean(se1^2,na.rm=T)*mean(se2^2,na.rm=T))
    r=cov_b1_b2/sqrt(var_b1*var_b2)
    if(jackknife){
    r_jack=c()
    n=length(b1)
    for(k in 1:n) {
       b1_jack=b1[-k];se1_jack=se1[-k];var_b1_jack=var(b1_jack,na.rm=T)-mean(se1_jack^2,na.rm=T)
       b2_jack=b2[-k];se2_jack=se2[-k];var_b2_jack=var(b2_jack,na.rm=T)-mean(se2_jack^2,na.rm=T)
 	if(var_b1_jack<0){
	var_b1_jack=var(b1_jack,na.rm=T)
}
if(var_b2_jack<0){
	var_b2_jack=var(b2_jack,na.rm=T)
}
       rp_jack=rp[-k];
       cov_e1_jack_e2_jack=mean(rp_jack,na.rm=T)*sqrt(mean(se1_jack^2,na.rm=T)*mean(se2_jack^2,na.rm=T))
       cov_b1_b2_jack=cov(b1_jack,b2_jack,use="complete.obs")-cov_e1_jack_e2_jack
       r_tmp=cov_b1_b2_jack/sqrt(var_b1_jack*var_b2_jack)
       r_jack=c(r_jack,r_tmp)
    }
    r_mean=mean(r_jack,na.rm=T)
    idx=which(is.na(r_jack))
    if(length(idx)>0){
        se_r=sqrt((n-1)/n*sum((r_jack[-idx]-r_mean)^2))
    }else{
        se_r=sqrt((n-1)/n*sum((r_jack-r_mean)^2))
    }

    res<-cbind(r,se_r)
    }else{
	    res=r
    }
    return(res)
}

calcu_eff_diff<-function(b1,se1,b2,se2,theta){
	d=b1-b2
	var=se1^2+se2^2-2*theta*se1*se2
	chisq=d^2/var
	p=pchisq(chisq,1,lower.tail=F)
	res<-cbind(d,var,chisq,p)
}

caucu_theta=function(beta1,se1,beta2,se2){
	var_e1=mean(se1^2);var_e2=mean(se2^2)
	cov_b1_b2_b1=cov(beta1-beta2,beta1)
	cov_b1_b2_b2=cov(beta1-beta2,beta2)
	rho=(1/2*(var_e1+var_e2+cov_b1_b2_b2-cov_b1_b2_b1))/sqrt(var_e1*var_e2)
	return(rho)
}

calcu_b_meta<-function(beta,SE,theta=diag(length(beta))){
	S=SE%*%t(SE)*theta

	Imat=matrix(1,nrow=dim(S)[1],ncol=1)
	eigenmin=min(eigen(S)$values)
	D=eigen(S)$values
	U=eigen(S)$vectors
	idx=which(abs(D)<1.0e-06)
	Dinv=vector(mode="numeric",length=length(D))
	if(length(idx)>0){
	   Dinv[idx]=1.0e-06
	   Dinv[-idx]=1/D[-idx]
	}else{
	    Dinv=1/D
	}
	InvD=diag(Dinv)
	Vinv=U%*%InvD%*%t(U)

	b_meta=solve(t(Imat)%*%Vinv%*%Imat)%*%t(Imat)%*%Vinv%*%beta
	se_meta=sqrt(solve(t(Imat)%*%Vinv%*%Imat))
        p_meta=pchisq((b_meta/se_meta)^2,1,lower.tail=F)
	return(c(b_meta,se_meta,p_meta))
}

Decoupled<-function(beta,SE,theta){
    omega=diag(SE)%*%theta%*%diag(SE)
    Imat=matrix(1,nrow=dim(omega)[1],ncol=1)
    eigenmin=min(eigen(omega)$values)
    D=eigen(omega)$values
    U=eigen(omega)$vectors
    idx=which(abs(D)<1.0e-06)
    Dinv=vector(mode="numeric",length=length(D))
    if(length(idx)>0){
       Dinv[idx]=1.0e-06
       Dinv[-idx]=1/D[-idx]
    }else{
        Dinv=1/D
    }
    InvD=diag(Dinv)
    Vinv=U%*%InvD%*%t(U)

    omega_decompled=solve(diag(as.numeric(t(Imat)%*%Vinv)))
    SE_decompled=diag(sqrt(omega_decompled))
    beta_decompled=beta
    res=calcu_b_meta(beta_decompled,SE_decompled,diag(nrow=length(beta_decompled)))
    return(res)
}

pi.0.est<-function(p){
    m <- length(p)
    n <- 2*m
    pa <- p
    pa[(m+1):n] <- 2-p
    h <- 0.9*min(sd(pa),(quantile(pa,.75)-quantile(pa,.25))/1.34)*n^(-1/5)
    dens <- n^-1*h^-1*(2*pi)^-(1/2)*sum(exp(-(1/2)*((1-pa)/h)^2))
    pi.0.hat <- 2*dens
    return(pi.0.hat)
}

create_z=function(y,x,f=7){
  sd=sd(y, na.rm=T)
  mean=mean(y,na.rm=T)
  y[which(abs(y-mean)>sd*f)]=NA;
  b=coefficients(summary(lm(y~., data=x)))
  e=(y-as.matrix(x)%*%b[-1,])[,1]
  return((e-mean(e,na.rm=T))/sd(e,na.rm=T))
}

invNorm=function(x){
  sd=sd(x, na.rm=T)
  mean=mean(x,na.rm=T)
  x[which(abs(x-mean)>sd*5)]=NA;
  return((x-mean(x,na.rm=T))/sd(x,na.rm=T))
}

RINT=function(x){
    return(qnorm((rank(x,na.last="keep",ties.method ="random")-3/8)/sum(!is.na(x))))
}


rm_mhc=function(smr,mhcStart=25000000,mhcEnd=36000000){
	idx=which(smr$ProbeChr==6 & ((smr$Probe_bp<=mhcEnd & smr$Probe_bp>=mhcStart) | (smr$topSNP_bp<=mhcEnd & smr$topSNP_bp>=mhcStart)))
	if(length(idx)>0){
	    smr=smr[-idx,]
	}
	return(smr)
}


rm_mhc_hg38=function(smr,mhcStart=28510120,mhcEnd=33480577){
    idx=which(smr$ProbeChr==6 & ((smr$Probe_bp<=mhcEnd & smr$Probe_bp>=mhcStart) | (smr$topSNP_bp<=mhcEnd & smr$topSNP_bp>=mhcStart)))
    if(length(idx)>0){
        smr=smr[-idx,]
    }
    return(smr)
}

rm_mhc_hg19=function(smr,mhcStart=28477797,mhcEnd=33448354){
    idx=which(smr$ProbeChr==6 & ((smr$Probe_bp<=mhcEnd & smr$Probe_bp>=mhcStart) | (smr$topSNP_bp<=mhcEnd & smr$topSNP_bp>=mhcStart)))
    if(length(idx)>0){
        smr=smr[-idx,]
    }
    return(smr)
}





get_com_SMR<-function(file1,file2){
	source("/shares/compbio/Group-Yang/t.qi/bin/CommonFunc.r")
	data1=fread(file1,head=T,stringsAsFactors=F,data.table=F)
	data2=fread(file2,head=T,stringsAsFactors=F,data.table=F)
	data1=rm_mhc(data1)
	data2=rm_mhc(data2)
	index=match(data1$probeID,data2$probeID,nomatch=0)
	data1_com=data1[which(index!=0),]
	data2_com=data2[index,]
	data_list=list(data1=data1_com,data2=data2_com)
	return(data_list)
}



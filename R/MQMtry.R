scan <- function(cross,cofactors,step.size=3,window.size=10){
	
	cross <- calc.genoprob(cross,step=step.size)
	cross <- est.rf(cross)
	num_markers <- 0
	step <- 0
	off_end <- 0
	markerdistances <- NULL
	geno <- NULL
	
	num_genotypes <- dim(cross$geno[[1]]$prob)[3]
	num_ind <- dim(cross$geno[[1]]$prob)[1]
	pheno <- cross$pheno[[pheno.col]]
	markerrf <- cross$rf
	
	for(x in 1:nchr(cross)){
		attr(cross$geno[[x]]$prob,"map")
		markerdistances <- c(markerdistances,as.vector(attr(cross$geno[[1]]$prob,"map")))
		num_markers = num_markers+length(attr(cross$geno[[x]]$prob,"map"))
		geno <- cbind(geno,cross$geno[[i]]$data)
	}

	probabilitymatrix <- array(rep(0,num_ind*num_markers*num_genotypes),c(num_ind,num_markers,num_genotypes))

	for(ind in 1:num_ind){
		mar <- 1	#index of marker inside the chromosome
		chr <- 1	#which chromosome we at
		for(mark in 1:num_markers){
			if(mark >  sum(nmar(cross)[1:chr])){
				#cat("Switching to next chr after marker",mark,"\n")
				mar <- 1
				chr <- chr+1
			}else{
				#cat(chr,ind,mark,geno,":",chr,ind,mar,geno,"\n")
			}
			for(geno in 1:num_genotypes){
				probabilitymatrix[ind,mark,geno] <- cross$geno[[chr]]$prob[ind,mar,geno]
			}
			mar = mar+1;
		}
	}

	#NOW WE KNOW ENOUGH IN THEORY TO EXECUTE MQM
	#-location of markers and their names and probabilities
	#-number of genotypes
	#-number of individuals
	
	.C("R_testScan",num_ind,num_markers,num_genotypes,trait,markerdistances,geno,probabilitymatrix,markerrf,cofactors,step.size,window.size)
	
}


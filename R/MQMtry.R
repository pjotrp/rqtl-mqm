scan <- function(cross,cofactors,step.size=3,window.size=10){

	library(qtl)
	data(hyper)
	step.size=3
	pheno.col=1
	window.size=10
	cross <- hyper
	cross <- fill.geno(cross)
	cross <- calc.genoprob(cross,step=step.size)
	cross <- est.rf(cross)
	num_markers <- 0
	step <- 0
	off_end <- 0
	markerdistances <- NULL
	genomatrix <- NULL
	
	num_genotypes <- dim(cross$geno[[1]]$prob)[3]
	num_ind <- nind(cross)
	trait <- cross$pheno[[pheno.col]]
	markerrf <- cross$rf
	cofactors <- rep(0,num_markers)
	markers <- NULL
	
	for(x in 1:nchr(cross)){
		attr(cross$geno[[x]]$prob,"map")
		markerdistances <- c(markerdistances,as.vector(attr(cross$geno[[x]]$prob,"map")))
		num_markers = num_markers+length(attr(cross$geno[[x]]$prob,"map"))
		markers <- c(markers,length(attr(cross$geno[[x]]$prob,"map")))
		genomatrix <- cbind(genomatrix,cross$geno[[x]]$data)
	}

	probabilitymatrix <- array(rep(0,num_ind*num_markers*num_genotypes),c(num_ind,num_markers,num_genotypes))

	for(ind in 1:num_ind){
		mar <- 1	#index of marker inside the chromosome
		chr <- 1	#which chromosome we at
		for(mark in 1:num_markers){
			if(mark >  sum(markers[1:chr])){
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
		#cat("Individual",ind,"transfered\n")
	}
	
	
	#NOW WE KNOW ENOUGH IN THEORY TO EXECUTE MQM
	#-location of markers and their names and probabilities
	#-number of genotypes
	#-number of individuals
	
	.C("R_testScan",num_ind,num_markers,num_genotypes,trait,markerdistances,genomatrix,probabilitymatrix,markerrf,cofactors,step.size,window.size)
	
}


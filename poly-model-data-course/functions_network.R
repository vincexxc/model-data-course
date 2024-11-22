corrected_laplacian<-function(network){
	lap<-laplacian_matrix(network)
	lap+diag(degree(network,mode="in"))
}

trophic_levels<-function(network) { #based on MacKay et al. 2020
	lap<-corrected_laplacian(network)
	imbalance<-degree(network,mode="in")-degree(network,mode="out")
	inv(as.matrix(lap)) %*% imbalance
}

as_Community<-function(network,net_name="."){#takes a directed network (with species names) and make it a Community object for cheddar
	names_net<-V(network)$name
	nodes<-data.frame("node"=names_net,row.names= names_net)
	properties=list("title"=net_name)
	tlinks<-as_long_data_frame(network)
	trophic.links<-as.matrix(tlinks[,c("from_name","to_name")])
	colnames(trophic.links)<-c("resource", "consumer")
	Community(nodes, properties, trophic.links)
}
#example: ShortestTrophicLevel(as_Community(network,"WTF"))

layout_as_food_web<-function(network){#adapted from Jon Borelli's code https://assemblingnetwork.wordpress.com/2013/04/03/more-food-web-plotting-with-r/
	l<-length(V(network))
	lay<-matrix(nrow=l,ncol=2) 
	lay[,1]<-layout_with_graphopt(network)[,1]
	lay[,2]<-TrophicLevels(as_Community(network))[,5]-1
	lay
}

layout_as_food_web2<-function(network){#adapted from Jon Borelli's code https://assemblingnetwork.wordpress.com/2013/04/03/more-food-web-plotting-with-r/
	l<-length(V(network))
	lay<-matrix(nrow=l,ncol=2) 
	lay[,1]<-layout_with_graphopt(network)[,1]
	lay[,2]<-trophic_levels(network)
	lay
}

layout_as_food_web3<-function(network){#adapted from Jon Borelli's code https://assemblingnetwork.wordpress.com/2013/04/03/more-food-web-plotting-with-r/
	l<-length(V(network))
	lay<-matrix(nrow=l,ncol=2) 
	lay[,1]<-layout_with_graphopt(network)[,1]
	lay[,2]<-alt_TL(network)[,3]
	lay
}

alt_TL<-function(network){ #yet another implementation of trophic levels
	undir_net<-as.undirected(network)
	basals<-which(degree(network,mode="in")==0)
	dist_mat<-t(distances(network,v=which(degree(network,mode="in")==0),mode="out"))
	s<-dim(dist_mat)[1]
	shortest_chain<-sapply(1:s,function(x) min_without_inf(dist_mat[x,]))+1
	longest_chain<-sapply(1:s,function(x) max_without_inf(dist_mat[x,]))+1
	average_chain<-sapply(1:s,function(x) mean_without_inf(dist_mat[x,]))+1
	if(!is.null(V(network)$name)){
		res<-data.frame("species"=V(network)$name,"shortest"=shortest_chain,"longest"=longest_chain,"average" = average_chain)
	}
	else{
		res<-data.frame("shortest"=shortest_chain,"longest"=longest_chain,"average" = average_chain)
	}
	res
}

min_without_inf<-function(vec){
	min(vec[!is.infinite(vec)])
}

max_without_inf<-function(vec){
	max(vec[!is.infinite(vec)])
}

mean_without_inf<-function(vec){
	mean(vec[!is.infinite(vec)])
}


cascade_matrix<-function(cc,nspecies){ #from Cohen-Newman-Briand's papers
	upper.tri(matrix(1,nrow=nspecies,ncol=nspecies), diag = FALSE)*matrix(rbinom(nspecies*nspecies,1,cc/nspecies),nrow=nspecies)
}

niche_matrix<-function(connectance,nspecies){ #Williams-Martinez model
	n_i<-runif(nspecies)
	r_i<-rbeta(nspecies,1,(1/(2*connectance))-1)
	r_i<-r_i*n_i
	c_i<-sapply(1:nspecies,function(x) runif(1,min=r_i[x]/2,max=n_i[x]))
	pred_function<-function(z_1,z_2){
		if((n_i[z_1]<=c_i[z_2]+0.5*r_i[z_2])&(n_i[z_1]>=c_i[z_2]-0.5*r_i[z_2])) {
			1
		}
		else{
			0
		}
	}
	mat<-sapply(1:nspecies,function(x) sapply(1:nspecies,function(y) pred_function(y,x)))
	list("matrix"=mat,"n"=n_i,"r"=r_i,"c"=c_i)
}


to_upper_triangular<-function(mat){
	upper.tri(mat, diag = FALSE)*mat
}

make_alluvial_2<-function(classif1,classif2,name1,name2){
	A <- as.data.frame(table(classif1,classif2))
	colnames(A) = c(name1,name2,"Freq")
	w   <- which(A$Freq != 0)
	A <- A[w,]
	alluvial(A[,c(1,2)],freq = A$Freq)
}

FW_interaction_from_predation<-function(mat,rho){
	n<-dim(mat)[1]
	fill<-rnorm_multi(n*(n-1)/2,2,r=rho)
	ut<-matrix(0,n,n)
	ut[lower.tri(ut,diag = FALSE)]<-fill[,1]
	ut<-mat*t(ut)
	lt<-matrix(0,n,n)
	lt[lower.tri(lt, diag = FALSE)]<-fill[,2]
	lt<-(t(mat))*lt
	ut+lt
}

spectral_clustering <- function(graph, nb_cluster, normalized = TRUE) {#from J. Chiquet's git page https://jchiquet.github.io/MAP566/docs/mixture-models/map566-lecture-graph-clustering-part1.html
  
  ## Compute Laplcian matrix
  L <- laplacian_matrix(graph, normalized = normalized) 
  ## Generates indices of last (smallest) K vectors
  selected <- rev(1:ncol(L))[1:nb_cluster] 
  ## Extract an normalized eigen-vectors
  U <- eigen(L)$vectors[, selected, drop = FALSE]  # spectral decomposition
  U <- sweep(U, 1, sqrt(rowSums(U^2)), '/')    
  ## Perform k-means
  res <- kmeans(U, nb_cluster, nstart = 40)$cl
  
  res
}

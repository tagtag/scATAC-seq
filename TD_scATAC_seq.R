TD_scATAC_seq <- function(x_files,L=10,nseed=0)
{
    require(Matrix)
    require(irlba)
    require(umap)
    cat("computing individual SVD\n")
    SVD_all <- NULL
    for (i in c(1:length(x_files)))
    {
        cat(i, " ")
        load(x_files[i])
        A <- t(t(x_all)/colSums(x_all))
        SVD <- irlba(A,L)
        a <- which(SVD$u!=0,arr.ind=T)
        SVDp <- spMatrix(dim(SVD$u)[1],dim(SVD$u)[2],i=a[,1],j=a[,2],x=SVD$u[a])
        SVD_all <- cbind(SVD_all,SVDp)
    }
    cat("\n computing total SVD\n")
    SVD_all_all <- irlba(SVD_all,L)
    a <- which(SVD_all_all$u!=0,arr.ind=T)
    SVD_all_all_p <- spMatrix(dim(SVD_all_all$u)[1],dim(SVD_all_all$u)[2],i=a[,1],j=a[,2],x=SVD_all_all$u[a])
    cat("computing projection\n")
    X<-NULL
    for (i in c(1:length(x_files)))
    {
    cat(i, " ")
    load(x_files[i])
    X <- rbind(X, t(x_all) %*% SVD_all_all_p)
    }
    cat("\n computing umap\n")
    set.seed(nseed)
    umap.defaults$n_neighbors <-30
    UMAP <- umap(X,config=umap.defaults)
    return(UMAP)
}
x_files<-c("../x_all_5","../x_all_6","../x_all_7","../x_all_8","../x_all_9","../x_all_10","../x_all_11","../x_all_12")


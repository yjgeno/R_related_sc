#experiment: row shared genes, col cells
object1 <- matrix(1:120, nrow = 10, ncol = 12)
object2 <- matrix(1:80, nrow = 10, ncol = 8)

library('irlba')
mat3 <- crossprod(x = object1, y = object2) #t(object1)%*%(object2)
cca.svd <- irlba(A = mat3, nv = min(nrow(mat3), ncol(mat3))-1)
cca.data <- rbind(cca.svd$u, cca.svd$v) # stack left and right singular vectors for cells

#flip certain cols 
cca.data <- apply(X = cca.data, MARGIN = 2, FUN = function(x) {
    if (sign(x[1]) == -1) {x <- x * -1}
    return(x)})

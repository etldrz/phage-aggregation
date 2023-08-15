source('single_patch.R')

one.wrong <- matrix(ncol=2, cbind(c(6, 4, 5), c(1, 5, 4)))

corrected <- swap(one.wrong)

all.equal(matrix(ncol=2, cbind(c(4, 6, 5), c(1, 5, 4))), corrected, check.attributes=F)

noPoolChanges = matrix(0, nrow=21, ncol=9)
# Gets rid of one G
noPoolChanges[c(6, 7, 9),1] = -1
# Adds one G
noPoolChanges[16,1] = 1
# Gets rid of one S
noPoolChanges[c(8, 10),2] = -1
# Adds one S
noPoolChanges[17,2] = 1
# Gets rid of one A
noPoolChanges[c(6, 8, 11),3] = -1
# Adds one A
noPoolChanges[c(1, 18),3] = 1
# Gets rid of one B
noPoolChanges[c(7, 12),4] = -1
# Adds on B
noPoolChanges[c(2, 19),4] = 1
# Gets rid of one G+A
noPoolChanges[c(3, 13),5] = -1
# Adds one G+A
noPoolChanges[6,5] = 1
# Gets rid of one G+B
noPoolChanges[c(4, 14),6] = -1
# Adds on G+B
noPoolChanges[7,6] = 1
# Gets rid of one S+A
noPoolChanges[c(5, 15),7] = -1
# Adds on S+A
noPoolChanges[8,7] = 1
# Removes one RA
noPoolChanges[c(1, 20),8] = -1
# Removes one RB
noPoolChanges[c(2, 21),9] = -1

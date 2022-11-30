
##############################################
# Additional
# Oct 2022
##############################################


# A) #########################################
# End of S4_TRaitDataIII.R -------------------
# I perfomed some permanovas, but it was unclear whether
# they are robust to unbalanced data and to differences
# in the dispersion of points.

# Sources:
# https://www.researchgate.net/post/Permanova_using_Adonis2_with_unbalanced_design
# https://sites.google.com/site/mb3gustame/warnings/warning-multiple-testing
# https://github.com/vegandevs/vegan/issues/344 (UNBALANCED)

##############################################
# PERMANOVA:----------------------------------
# Used to compare groups of objects and test the null hypothesis 
# that the centroids and dispersion of the groups as defined 
# by measure space are equivalent for all groups.

identical(colnames(dist_mat1), rownames(tt)) #TRUE
grouping <- as.factor(tt$Status)

pm1 <- adonis(dist_mat1~grouping) #adonis2?
pm1 # R2 = 0.089
#pm2 <- adonis(coords~grouping, method = "euclidean")
# R-squared of 0.089*100= 8.9%

grouping2 <- as.factor(tt$SourceII)

pm2 <- adonis(dist_mat1~grouping2) # adonis 2?
pm2 # R2 0.17439
# R-squared of 0.17439*100= 17.43%

pairwise.adonis(coords, grouping, sim.function='vegdist',sim.method='euclidian')
pairwise.adonis(coords, grouping2, sim.function='vegdist',sim.method='euclidian')
# OK, so overall we see that the group Aquaculture&Sportfishing is the most
# different in composition. (But sample size)
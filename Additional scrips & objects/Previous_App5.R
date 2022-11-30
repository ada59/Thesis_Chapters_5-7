

# Test observer (fishtrait_fill and fishmorph)
# Need to load and subset FMorph first as in the equivalent appendix
# Previous:---------------------------------------------------------------------------------
load("fishtrait_fill.RData")

fishtrait_fill$Genus <- str_split_fixed(fishtrait_fill$Genus, " ", 2)[,1]
fishtrait_fill$Observer <- "AFE"
fishtrait_fill <- hablar::retype(fishtrait_fill)
str(fishtrait_fill)

FMorph_AFEII <- bind_rows(F_Morph_subset, fishtrait_fill) # combine both datasets
picture_traits_combined <- FMorph_AFEII[,c(7:15)]
rownames(picture_traits_combined) <- FMorph_AFEII$Genus.species
picture_traits_combined_II <- missForest(as.matrix(picture_traits_combined))

PCAFMAFEII <- prcomp(picture_traits_combined_II$ximp, scale=T) # PCA combined

(screep <- fviz_eig(PCAFMAFEII, addlabels=TRUE, hjust = -0.3))

(bI_PCAFMAFEII <- fviz_pca_ind(PCAFMAFEII,
                               title = "",
                               label="none",
                               habillage = FMorph_AFEII$Genus, # color by groups
                               addEllipses = TRUE, # TRUE for concentration ellipses
                               legend.title="Genus",
                               ellipse.type = "convex",
                               alpha.ind = 0.5))
observer_groups <- FMorph_AFEII$Observer
observer_groups[observer_groups=="FM"] <- ""
(bI_PCAFMAFEII <- bI_PCAFMAFEII + geom_text(aes(label = observer_groups),
                                            alpha = 0.5, size = 2, nudge_y = 0.1, show.legend = FALSE))

(bI_PCAFMAFEII_3_4 <- fviz_pca_ind(PCAFMAFEII,
                                   title = "",
                                   axes = c(3,4),
                                   label="none",
                                   habillage = FMorph_AFEII$Genus, # color by groups
                                   addEllipses = TRUE, # TRUE for concentration ellipses
                                   legend.title="Genus",
                                   ellipse.type = "convex",
                                   alpha.ind = 0.5))
observer_groups <- FMorph_AFEII$Observer
observer_groups[observer_groups=="FM"] <- ""
(bI_PCAFMAFEII_3_4 <- bI_PCAFMAFEII_3_4 + geom_text(aes(label = observer_groups),
                                                    alpha = 0.5, size = 2, nudge_y = 0.1, show.legend = FALSE))


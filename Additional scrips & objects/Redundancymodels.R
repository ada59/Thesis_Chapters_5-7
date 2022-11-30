##########################################################################################
# Observed trends (R0 vs T0):-------------------------------------------------------------
# Historical conservative:----------------------------------------------------------------
hnc_1R0 <- lm(R0 ~ T0, data=hnc)
summary(hnc_1R0)
hnc_2R0 <- lm(R0 ~ log(T0+1), data=hnc)

AIC(hnc_1R0, hnc_2R0) # Linear

# Historical broad:-----------------------------------------------------------------------
hnb_1R0 <- lm(R0 ~ T0, data=hnb)
summary(hnb_1R0)
hnb_2R0 <- lm(R0 ~ log(T0+1), data=hnb)

AIC(hnb_1R0, hnb_2R0) # Linear

# Current native:-------------------------------------------------------------------------
contN_1R0 <- lm(R0 ~ T0, data=contN)
summary(contN_1R0)
contN_2R0 <- lm(R0 ~ log(T0+1), data=contN)

AIC(contN_1R0, contN_2R0) # Linear

# Current All:----------------------------------------------------------------------------
contAll_1R0 <- lm(R0 ~ T0, data=contAll)
summary(contAll_1R0)
contAll_2R0 <- lm(R0 ~ log(T0+1), data=contAll)

AIC(contAll_1R0, contAll_2R0) # Linear
# NOTE: A linear regression is a better fit in all cases ( R0 vs T0)

##########################################################################################
# Assumptions: ---------------------------------------------------------------------------
# Histograms:-----------------------------------------------------------------------------
par(mfrow = c(2, 2))
hist(hnc_1R0$residuals)
hist(hnb_1R0$residuals)
hist(contN_1R0$residuals)
hist(contAll_1R0$residuals)
dev.off()

# Plots:----------------------------------------------------------------------------------
par(mfrow = c(2, 2), mar = c(2, 2, 2, 2))
plot(hnc_1R0)
plot(hnb_1R0)
plot(contN_1R0)
plot(contAll_1R0)
dev.off()

##########################################################################################
# Plot trends:----------------------------------------------------------------------------

(p1R0 <- ggplot(div, aes(x=T0, y=R0, color=Period))+
   geom_point(aes(color=Period))+
   geom_smooth(method=lm, formula = y ~ x)+
   theme_classic()+
   scale_color_manual(values=colorBlind4)+
   facet_wrap(~Period)) 

ggsave(p1R0, file= paste0(plot_dir, "/SM/Thesis/RedundancyLinear4.jpg"), width = 8, height = 6)



(nullR0 <- ggplot(data=F0null, aes(x=as.factor(T0_n), y=R0_n))+
    geom_boxplot(alpha=0.5, color="gray49")+
    stat_boxplot(geom ='errorbar', color="gray49") + 
    xlab("SR")+
    ylab("R0")+
    theme_minimal()) # null SR vs R0



(p1nullR0 <- nullR0 + geom_point(data = div, 
                                 aes(x=T0, y=R0, color=Period), size=2, alpha=0.5)+
    geom_smooth(data = div, aes(x=T0, y=R0, color=Period), method=lm, formula = y ~ x, se=TRUE)+
    scale_color_manual(values=colorBlind4)+
    facet_wrap(~Period)) # SM plot, trait redundancy
(p1nullR0_3 <- nullR0 + geom_point(data = div_subset_main, 
                                   aes(x=T0, y=R0, color=Period), size=2, alpha=0.5)+
    geom_smooth(data = div_subset_main, aes(x=T0, y=R0, color=Period), method=lm, formula = y ~ x, se=TRUE)+
    scale_color_manual(values=colorBlind3)+
    theme(legend.position = "none",
          axis.text.x = element_text(size=8))+
    facet_wrap(~Period)) # Main ms plot, trait redundancy


ggsave(p1nullR0, file= paste0(plot_dir, "/Null-ObservedR0.jpg"), width = 8, height = 6) 

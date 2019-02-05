library(ggplot2)
library(scales)
nproc = 2**seq(5,10)
meansSeq = c(2800,	1400,	908.9)
errorSeq = c(1800,	1000,	458.034)
fillSeq = rep("seq", times=3)
meansInv = c(204.13,	114.9,	56.3,	45.1,	29.6,	25.3)
errorInv = c(113.65,	37.27,	12.19,	15.31, 4.22,	4.49)
fillInv = rep("inv", times=6)
meansEpsBE = c(235.3,	76.4,	73.7,	44.5,	29.8,	22.5)
errorEpsBE = c(140.89, 21.98, 27.50, 10.66, 4.89,	4.09)
fillEpsBE = rep("be", times=6)
meansEpsRMSD = c(111.8,	84.5,	52.9,	42.1,	31.9,	24.5)
errorEpsRMSD = c(51.66,	29.69, 20.27,	11.22,	3.47,	2.95)
fillEpsRMSD = rep("rmsd", times=6)
errorCorrected = meansInv-errorInv
errorCorrected[1] = 1
dd <- data.frame(nproc, meansSeq, meansInv, meansEpsBE, meansEpsRMSD, errorCorrected, errorSeq, errorInv,errorEpsBE, fillSeq, fillInv, fillEpsRMSD, fillEpsBE)

ggplot(dd) +
  geom_ribbon(aes(ymax=meansSeq+errorSeq, ymin = meansSeq-errorSeq, x=nproc),fill="red", alpha=0.4)+
  geom_ribbon(aes(ymax=meansInv+errorInv, ymin = errorCorrected, x=nproc), fill = "blue", alpha=0.2)+
  geom_ribbon(aes(ymax=meansEpsRMSD+errorEpsRMSD, ymin = meansEpsRMSD-errorEpsRMSD, x=nproc), fill = "orange", alpha=0.5)+
  geom_ribbon(aes(ymax=meansEpsBE+errorEpsBE, ymin = meansEpsBE-errorEpsBE, x=nproc), fill = "green", alpha=0.3)+
  geom_point(aes(nproc, meansSeq, colour=fillSeq), size = 3, show.legend = TRUE ) +
  geom_point(aes(nproc, meansInv, colour=fillInv), size = 3, show.legend = TRUE ) +
  geom_point(aes(nproc, meansEpsBE, colour=fillEpsBE), size = 3, show.legend = TRUE) +
  geom_point(aes(nproc, meansEpsRMSD, colour=fillEpsRMSD), size = 3, show.legend = TRUE ) +
  # coord_trans(x="log2", y="log10")
  # show exponents
  scale_x_continuous(trans = log2_trans())+
                     # breaks = trans_breaks("log2", function(x) 2^x), 
                     # labels = trans_format("log2", math_format(2^.x))) +
  scale_y_continuous(trans = log10_trans())+
                     # breaks = trans_breaks("log10", function(x) 10^x),
                     # labels = trans_format("log10", math_format(10^.x)))
  guides(colour=guide_legend(title="Method"))+
  guides(fill=FALSE)+
  #guides(alpha=FALSE)+
  # guides(size=FALSE)+
  scale_colour_manual(name="Method",
                      breaks=c("seq", "inv", "be", "rmsd"),
                      values = c("green", "blue", "orange", "red"),
                      labels=c("Sequential", "Inversely", "Epsilon BE", "Epsilon RMSD"))+
  labs(title = "First binding event",
       y = "Steps", x = "Number of processors") +
  theme(title = element_text(size=16))+
  theme(legend.text = element_text(size=14)) +
  theme(axis.text = element_text(size=12))+
  theme(legend.title = element_text(face="bold"))+
  theme(legend.background = element_rect(colour = "black")) 
ggsave("plot4DAJ.png",device="png", width = 20, height = 10, units = "cm", dpi=300)

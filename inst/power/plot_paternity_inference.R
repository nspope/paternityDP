library(ggplot2)
library(dplyr)
library(cowplot)

nicely_format <- function(x, signif=2)
{
  x <- ifelse(x > 0, as.character(round(x, 2)), NA)
  x[x == '0'] <- '< 0.01'
  x
}

load("test_paternity_inference.realistic.RData")

power %>% dplyr::group_by(n_father, truth, mistyping) %>% 
  dplyr::summarise(
    post_prob=mean(post_prob), 
    color=ifelse(unique(truth)==unique(n_father),"light","dark")
  ) %>% 
  dplyr::mutate(
    error_rates = paste0("True error rates = ", mistyping)
  ) %>%
  ggplot(aes(x=truth, y=n_father)) + 
  geom_tile(aes(fill=post_prob)) +
  geom_text(aes(label=ifelse(truth==n_father, nicely_format(post_prob), NA), color=color), size=3, family="Garamond") +
  scale_fill_gradient("Average\nPost. Prob.", low="white", high="navyblue") + 
  scale_color_manual(values=c("light"="white","dark"="black")) +
  guides(color='none') +
  scale_y_continuous(expand=c(0,0), breaks=seq(1:8), limits=c(0.5,6.5)) +
  scale_x_continuous(expand=c(0,0), breaks=seq(1:6)) +
  theme_cowplot() + 
  theme(
    #legend.position=c(0.05,0.7),
    #legend.title=element_blank(),
    #plot.title=element_text(face="bold", hjust=1),
    plot.background=element_rect(fill="white"), 
    text=element_text(family="Garamond")
  ) +
  facet_grid(~error_rates) +
  #ggtitle("A.") +
  xlab("True # fathers") + ylab("Est. # fathers") -> plt_A

power %>% dplyr::group_by(truth, seed, mistyping) %>% 
  dplyr::summarise(
    post_dropout=mean(post_dropout), 
    post_mistyping=mean(post_mistyping) 
  ) %>%
  dplyr::mutate(
    error_rates = paste0("True error rates = ", mistyping)
  ) %>%
  ggplot(aes(x=truth, group=truth)) + 
  geom_hline(data=data.frame(error_rates = paste0("True error rates = ", c(0.01, 0.05, 0.10)),
                             mistyping = c(0.01, 0.05, 0.10)),
             aes(yintercept=mistyping), linetype=1, size=0.5) +
  geom_boxplot(aes(x=truth-0.18, y=post_mistyping, fill="Mistyping"), width=0.3, outlier.size=0.75) + #, outlier.shape=NA) +
  geom_boxplot(aes(x=truth+0.18, y=post_dropout, fill="Dropout"), width=0.3, outlier.size=0.75) + #, outlier.shape=NA) +
#  geom_violin(aes(x=truth-0.13, y=post_mistyping, fill="Mistyping error"), width=0.2) +
#  geom_violin(aes(x=truth+0.13, y=post_dropout, fill="Dropout error"), width=0.2) +
  scale_y_continuous(breaks=seq(0, 0.2, 0.1), labels=c(" 0.0", " 0.1", " 0.2")) +
  scale_x_continuous(breaks=seq(1:6)) +
  scale_fill_manual("Error class", values=c("Mistyping"="dodgerblue", "Dropout"="coral")) + 
  theme_cowplot() + 
  theme(
    #legend.position=c(0.05,0.7),
    #legend.title=element_blank(),
    #plot.title=element_text(face="bold", hjust=1),
    plot.background=element_rect(fill="white"), 
#    axis.text=element_text(size=10),
    text=element_text(family="Garamond")
  ) +
  facet_grid(~error_rates) +
#  ggtitle("B.") +
  xlab("True # fathers") + ylab("Est. error rate") -> plt_B

png("test_paternity_inference.png", height=5.5, width=9, units="in", res=300)
cowplot::plot_grid(plt_A, plt_B, nrow=2, labels="AUTO", label_fontfamily="Garamond", label_size=16, align="v")
dev.off()

#ggsave("test_paternity_inference.png", plt_A, height=3, width=5, units="in", dpi=300)

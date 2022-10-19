mainpaper_plot
================
Xinran Wang
2022-10-19

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
    ## ✔ tibble  3.1.7     ✔ dplyr   1.0.9
    ## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
    ## ✔ readr   2.1.2     ✔ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library(ggplot2)
library(cowplot)
library(rlang)
```

    ## 
    ## Attaching package: 'rlang'

    ## The following objects are masked from 'package:purrr':
    ## 
    ##     %@%, as_function, flatten, flatten_chr, flatten_dbl, flatten_int,
    ##     flatten_lgl, flatten_raw, invoke, splice

Canonical Parameters GWAS simulation mode: fast, linear model: enet, SNP
model: 1pct, eQTL sample size: 250, GWAS sample size: 200000, h2ge:
0.0001, h2g: 0.1.

**Dataframe 1**

``` r
setwd("~/Local Documents/Software_Output/R_Studio_Output/TWAS/Results0818")
df_mem <- read_tsv("mem.tsv")
```

    ## Rows: 480 Columns: 9
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): gwas, snp_model, linear_model
    ## dbl (6): id, ngwas, MaxRSS, h2ge, h2g, nge
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
df_barchart_mem <- df_mem %>%
  group_by(ngwas, gwas) %>%
  summarise(ngwas = mean(ngwas), 
            n=n(), 
            mean.MaxRSS = mean(MaxRSS),
            sd = sd(MaxRSS),
            se = sd / sqrt(n))
```

    ## `summarise()` has grouped output by 'ngwas'. You can override using the
    ## `.groups` argument.

``` r
df_barchart_mem$ngwas <- factor(df_barchart_mem$ngwas, levels = c(50000, 100000, 500000), labels = c("FiftyK", "OneHundK", "FiveHundK"))
df_barchart_mem$gwas <- gsub('Fast mode', 'fast', df_barchart_mem$gwas)
df_barchart_mem$gwas <- gsub('Standard mode', 'std', df_barchart_mem$gwas)
```

**Dataframe 2**

``` r
setwd("~/Local Documents/Software_Output/R_Studio_Output/TWAS/Results0914")
df_summary <- read_tsv("all.summary.tsv")
```

    ## Rows: 115200 Columns: 20
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (3): gwas.sim, linear_model, snp_model
    ## dbl (17): sim, id, real.time, cpu.time, h2ge, nsnps, ngwas, nqtl, h2g, h2g.h...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
df_barchart_cputime <- df_summary %>%
  group_by(ngwas, gwas.sim) %>%
  summarise(ngwas = mean(ngwas), 
            n=n(), 
            mean.cpu.time = mean(cpu.time),
            sd = sd(cpu.time),
            se = sd / sqrt(n))
```

    ## `summarise()` has grouped output by 'ngwas'. You can override using the
    ## `.groups` argument.

``` r
df_qq <- df_summary %>% filter(gwas.sim == "fast", ngwas == 200000 , nqtl == 250, snp_model == "1pct", h2ge == 0) 
df_h2ge_power <- df_summary %>% filter(gwas.sim == "fast", ngwas == 200000, nqtl == 250, snp_model == "1pct")
df_barchart_cputime$ngwas <- factor(df_barchart_cputime$ngwas, levels = c(2500, 5000, 7000),  labels = c("2500", "5000", "7000"))
```

**Plot Settings** Color Preset

``` r
lm.group.colors <- c(enet="#F8766D", lasso="#7CAE00", ridge="#00BFC4", trueqtl="#C77CFF")
gwas.mode.group.colors <- c(fast="darkblue", std="steelblue3")
```

Param Preset

``` r
thold <- 0.05/22000
```

Plot Preset

``` r
font.size <- 10
text.angel <- 45
text.vjust <- 0.75
text.hjust <- 0.50
cowplot.margin <- margin(20, 0, 40, 0)
```

Legend Names

``` r
ngwas_names1 <- c('FiftyK', 'OneHundK', 'FiveHundK')
ngwas_names2 <- c("50K", "100K", "500K")
lm_names1 <- c('enet', 'lasso', 'ridge', 'trueqtl')
lm_names2 <- c("Elastic Net", "LASSO", "GBLUP", "True eQTL")
gwas_mode_names1 <- c('fast', 'std')
gwas_mode_names2 <- c("Fast", "Standard")
```

**Cow Plot**

``` r
plot_1 <- ggplot(df_qq, aes(sample = twas.z^2, color= linear_model, fill= linear_model, group= linear_model)) +
  geom_qq(distribution = stats::qchisq, dparams=list(df=1)) +
  geom_abline(slope=1, intercept=0, color="black") +
  labs(x = "Expected",
       y = "Observed") +
  theme(axis.text = element_text(size = font.size),
        legend.direction = "horizontal",
        legend.title=element_text(size=font.size),
        legend.text=element_text(size=font.size), 
        axis.text.x = element_text(angle =text.angel, vjust=text.vjust, hjust=text.hjust), 
        panel.background = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  guides(fill=guide_legend(title="Linear Model"),
         color=guide_legend(title="Linear Model", nrow = 1, byrow = TRUE)) + 
  scale_color_manual(values=lm.group.colors, breaks = lm_names1, labels = lm_names2) + 
  scale_fill_manual(values=lm.group.colors, breaks = lm_names1, labels = lm_names2) 

plot_2 <- ggplot(df_h2ge_power %>% mutate(power = as.numeric(twas.p < thold)), aes(x=as.factor(h2ge))) + 
  stat_summary(aes(y=power, group=linear_model, color=linear_model), fun.data = "mean_cl_boot", geom="point") + 
  stat_summary(aes(y=power, group=linear_model, color=linear_model), fun.data = "mean_cl_boot", geom="errorbar", fun.args = list(conf.int = 0.95), width = 0.2) + 
  stat_summary(aes(y=power, group=linear_model, color=linear_model), fun.data = "mean_cl_boot", geom="line") +
  ylim(0,1) +
  labs(x = bquote(""~~ h[GE]^2 ~~""), y = "Power") +
  theme(axis.text = element_text(size = font.size),
        legend.direction = "horizontal",
        legend.title=element_text(size=font.size),
        legend.text=element_text(size=font.size), 
        axis.text.x = element_text(angle =text.angel, vjust=text.vjust, hjust=text.hjust), 
        panel.background = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  guides(fill=guide_legend(title="Linear Model"),
         color=guide_legend(title="Linear Model", nrow = 1, byrow = TRUE)) + 
  scale_color_manual(values=lm.group.colors, breaks = lm_names1, labels = lm_names2) + 
  scale_fill_manual(values=lm.group.colors, breaks = lm_names1, labels = lm_names2) 

plot_3 <- ggplot(df_barchart_mem) +
  geom_bar(aes(x=ngwas, y=mean.MaxRSS/1000000, color = gwas, fill = gwas, group = gwas), 
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(x=ngwas, ymin=mean.MaxRSS-sd, ymax=mean.MaxRSS+sd, group = gwas), 
                colour="orange", position = position_dodge(width = 0.9), width=0.2, size=0.5) +
  labs(x = bquote("GWAS Sample Size"),
       y = bquote("Memory (MB)")) +
  ylim(0,9000) +
  scale_x_discrete(breaks=ngwas_names1, labels=ngwas_names2) +
  theme(axis.text = element_text(size = font.size), 
        axis.text.x = element_text(angle = 0, vjust=text.vjust, hjust=text.hjust), 
        legend.direction = "horizontal", 
        legend.title=element_text(size=font.size), 
        legend.text=element_text(size=font.size), 
        panel.background  = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  guides(fill=guide_legend(title="GWAS Mode"),
         color=guide_legend(title="GWAS Mode", nrow = 1, byrow = TRUE)) +
  scale_color_manual(values=gwas.mode.group.colors, breaks=gwas_mode_names1, labels=gwas_mode_names2) + 
  scale_fill_manual(values=gwas.mode.group.colors, breaks=gwas_mode_names1, labels=gwas_mode_names2)

prow <- plot_grid(
  plot_1 + theme(plot.margin = margin(6, 0, 40, 0), legend.position="none"),
  plot_2 + theme(plot.margin = margin(6, 0, 40, 0), legend.position="none"),
  plot_3 + theme(plot.margin = margin(6, 0, 40, 0), legend.position="none"),#legend.position=c(0.4, 0.8)),
  align = 'vh',
  labels = c("A", "B", "C"),
  hjust = -1,
  nrow = 1
)
prow
```

![](mainpaper_plot_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
legend1 <- get_legend(plot_1 + theme(legend.position = c(0.72, 2.4)))
legend3 <- get_legend(plot_3 +theme(legend.position = c(-0.6, 0.9)))
combineLegend <- plot_grid(legend1, legend3)

cow_plot <- plot_grid(prow, combineLegend, ncol = 1, rel_heights = c(1, .1))
# save_plot(file = "mainpaper_plot.pdf", plot = cow_plot, base_width = 17.4, base_height=6.5, units = "cm")
```

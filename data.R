library (tidyverse)
group_comparisons <- read_csv("group comparisons.csv")
df<- group_comparisons


##
# tt comes from `analyse.r`
## ind <- as.data.frame(tt)$Family=="Muribaculaceae"
## Muri_genus <- unique(as.data.frame(tt)$Genus[ind])
## Muri_genus <- Muri_genus[!is.na(Muri_genus)]
## is_ambiguous <- sapply(Muri_genus, function(g) {
##   g_ind <- as.data.frame(tt)$Genus==g
##   g_ind <- !is.na(g_ind) & g_ind
##   fam <- unique(as.data.frame(tt)$Family[g_ind])
##   length(fam)>1
## }
## )
## Muri_genus <- Muri_genus[!is_ambiguous]

## These are the genus' that Muribaculaceae
highlight <- list(Muribaculaceae = c(
  "Muribaculaceae uncultured bacterium", "Muribaculaceae Ambiguous_taxa",
  "uncultured Bacteroidales bacterium", "Muribaculum",
  "CAG-873", "mouse gut metagenome",
  "Gram-negative bacterium cTPY-13"))

library(ggrepel)
volcano <- function(de, highlight=list(),
            alpha=0.05, lfc=0.6,
            colmap = c(UP="forestgreen", DOWN="steelblue3"),
            title
            ) {
  if (missing(title)) {
    title <- deparse(substitute(de))
  }
  de$padj<- as.numeric(de$padj)
  de$diffexpressed <- "NO"
  de$diffexpressed[de$lfc > lfc & de$padj < alpha] <- "UP"
  de$diffexpressed[de$lfc < -lfc & de$padj < alpha] <- "DOWN"
  de$delabel <- NA
  de$delabel[de$diffexpressed != "NO"] <- de$level[de$diffexpressed != "NO"]
  de$shape <- factor("No", levels=c("No", names(highlight)))
  for (i in names(highlight)) {
    de$shape[de$level %in% highlight[[i]]] <- i
  }
  de$delabel[de$shape!="No"] <- de$level[de$shape!="No"]
  pl <- ggplot(data=de, aes(x=lfc, y=-log10(padj), col=diffexpressed, label=delabel, shape=shape)) +
    geom_point(size=3) + 
    theme_minimal() +
    geom_text_repel(show.legend=FALSE) +
    scale_color_manual(values=colmap, na.value="grey") +
    geom_vline(xintercept=c(-lfc, lfc), col="darkred", linetype="dotted") +
    geom_hline(yintercept=-log10(alpha), col="darkred", linetype="dotted") +
    guides(shape=ifelse(length(highlight)>0, "legend", "none")) +
    labs(color="Significant direction",
         shape="Category of interest",
         x="Log Fold Change",
         y="-Log10 adjusted p-value",
         title=title)
  print(pl)
  pl
}


df_glu_week2 <-df %>% filter( Week == 2 & Treatment == 'glucose')

df_glu_week2_phylum <- df_glu_week2 %>% filter( rank == "Phylum")

df_glu_week2_class <- df_glu_week2 %>% filter( rank == "Class")

df_glu_week2_order <- df_glu_week2 %>% filter( rank == "Order")

df_glu_week2_genus <- df_glu_week2 %>% filter( rank == "Genus")

###-----glucose Genus week 2####
de <- df_glu_week2_genus
volcano(df_glu_week2_genus, highlight=highlight)
# or volcano(df_glu_week2_genus, highlight=highlight, title="Glucose Week 12, Genus")
# omit 'highlight=...' to remove Muribaculaceae annotation
# can add lfc= and alpha= to change the default thresholds

de$padj<- as.numeric(de$padj)

ggplot(data=de, aes(x=lfc, y=padj)) + geom_point()
p <- ggplot(data=de, aes(x=lfc, y= -log10(padj))) + geom_point()
p
p <- ggplot(data=de, aes(x=lfc, y= -log10 (padj))) + geom_point() + theme_minimal()
p
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="grey") +
  geom_hline(yintercept=-log10(0.05), col="red")
p2
de$diffexpressed <- "NO"

de$diffexpressed[de$lfc > 0.6 & de$padj < 0.05] <- "UP"
de$diffexpressed[de$lfc < -0.6 & de$padj < 0.05] <- "DOWN"

p <- ggplot(data=de, aes(x=de$lfc, y=-log10(de$padj), col=diffexpressed)) + geom_point() + theme_minimal()
p

p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
p2

mycolors <- c("grey", "dodgerblue", "forestgreen")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)





###---- glucose genus week 12####

df_glu_week12 <-df %>% filter( Week == 12 & Treatment == 'glucose')

df_glu_week12_genus <- df_glu_week12 %>% filter( rank == "Genus")

volcano(df_glu_week12_genus, highlight=highlight)

de <- df_glu_week12_genus
de$padj<- as.numeric(de$padj)
str(de)

de$diffexpressed <- "NO"
View(de)
de$diffexpressed[de$lfc > 0.6 & de$padj < 0.05] <- "UP"
de$diffexpressed[de$lfc < -0.6 & de$padj < 0.05] <- "DOWN"

p <- ggplot(data=de, aes(x=de$lfc, y=-log10(de$padj)), col=diffexpressed) + geom_point() + theme_minimal()
p

p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
p2
mycolors <- c("grey", "dodgerblue", "forestgreen")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)

de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- de$level[de$diffexpressed != "NO"]


ggplot(data=de, aes(x=lfc, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("steelblue3", "grey", "darkred")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black")

###----  highsucralose genus week 12 ####

df_hs_week12 <-df %>% filter( Week == 12 & Treatment == 'high_sucralose')

df_hs_week12_genus <- df_hs_week12 %>% filter( rank == "Genus")
volcano(df_hs_week12_genus, highlight=highlight)

de <- df_hs_week12_genus
de$padj<- as.numeric(de$padj)
str(de)

de$diffexpressed <- "NO"
View(de)
de$diffexpressed[de$lfc > 0.6 & de$padj < 0.05] <- "UP"
de$diffexpressed[de$lfc < -0.6 & de$padj < 0.05] <- "DOWN"

p <- ggplot(data=de, aes(x=de$lfc, y=-log10(de$padj)), col=diffexpressed) + geom_point() + theme_minimal()
p

p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
p2

de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- de$level[de$diffexpressed != "NO"]


ggplot(data=de, aes(x=lfc, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("steelblue3", "grey", "darkred")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black")


###----  highsucralose genus week 2 ####

df_hs_week2 <-df %>% filter( Week == 2 & Treatment == 'high_sucralose')

df_hs_week2_genus <- df_hs_week2 %>% filter( rank == "Genus")
volcano(df_hs_week2_genus, highlight=highlight)

de <- df_hs_week2_genus
de$padj<- as.numeric(de$padj)
str(de)

de$diffexpressed <- "NO"
de$diffexpressed[de$lfc > 0.6 & de$padj < 0.05] <- "UP"
de$diffexpressed[de$lfc < -0.6 & de$padj < 0.05] <- "DOWN"

p <- ggplot(data=de, aes(x=de$lfc, y=-log10(de$padj), col=diffexpressed)) + geom_point() + theme_minimal()
p

p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
p2

de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- de$level[de$diffexpressed != "NO"]


ggplot(data=de, aes(x=lfc, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("steelblue3", "grey", "darkred")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black")


###----  low_sucralose genus week 2 ####

df_ls_week2 <-df %>% filter( Week == 2 & Treatment == 'low_sucralose')

df_ls_week2_genus <- df_ls_week2 %>% filter( rank == "Genus")

volcano(df_ls_week2_genus, highlight=highlight)

de <- df_ls_week2_genus
de$padj<- as.numeric(de$padj)
str(de)

de$diffexpressed <- "NO"
View(de)
de$diffexpressed[de$lfc > 0.6 & de$padj < 0.05] <- "UP"
de$diffexpressed[de$lfc < -0.6 & de$padj < 0.05] <- "DOWN"

p <- ggplot(data=de, aes(x=de$lfc, y=-log10(de$padj)), col=diffexpressed)) + geom_point() + theme_minimal()
p

p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
p2

de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- de$level[de$diffexpressed != "NO"]


ggplot(data=de, aes(x=lfc, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("steelblue3", "grey", "darkred")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black")


###----  low_sucralose genus week 12 ####

df_ls_week12 <-df %>% filter( Week == 12 & Treatment == 'low_sucralose')

df_ls_week12_genus <- df_ls_week12 %>% filter( rank == "Genus")
volcano(df_ls_week12_genus, highlight=highlight)

de <- df_ls_week12_genus
de$padj<- as.numeric(de$padj)
str(de)

de$diffexpressed <- "NO"
View(de)
de$diffexpressed[de$lfc > 0.6 & de$padj < 0.05] <- "UP"
de$diffexpressed[de$lfc < -0.6 & de$padj < 0.05] <- "DOWN"

p <- ggplot(data=de, aes(x=de$lfc, y=-log10(de$padj)), col=diffexpressed)) + geom_point() + theme_minimal()
p

p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
p2

de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- de$level[de$diffexpressed != "NO"]


ggplot(data=de, aes(x=lfc, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("steelblue3", "grey", "darkred")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black")



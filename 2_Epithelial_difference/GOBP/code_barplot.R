library(ggplot2)
df = read.table("./GOBP_GSEA_barplot.txt", header=T, sep="\t")
df$NES = round(df$NES, 4)
ggplot(df, aes(x=NES, y=reorder(Description, NES), fill=pvalue))+
geom_bar(stat='identity', width=0.6)+
scale_fill_gradient(low="#EE82EE", high="#FFD39B")+
theme_minimal()+geom_text(aes(label=NES), hjust=0)+
theme(text=element_text(family="Times"), axis.text.x=element_text(face="bold", size=14), axis.text.y=element_text(face="bold", size=14), axis.title.y=element_blank(), axis.title.x=element_text(face="bold", size=15))




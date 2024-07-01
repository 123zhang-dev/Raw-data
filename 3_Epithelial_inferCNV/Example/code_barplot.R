library(ggplot2)

df = read.table("./BP_HPV_neg_chr17p_del_barplot.txt", header=T, sep="\t")
ggplot(df, aes(x=Count, y=reorder(Term, Count), fill=PValue))+
geom_bar(stat="identity", width=0.5)+
scale_fill_gradient(low="dodgerblue", high="lightgrey")+
theme_bw()+
theme(axis.text.x=element_text(face="bold", size=12), axis.text.y=element_text(face="bold", size=12), axis.title.x=element_text(face="bold", size=15), axis.title.y=element_blank())






library(ggplot2)

df = read.table("./BP_Epithelial_subcluster_3_barplot.txt", header=T, sep="\t")
ggplot(df, aes(x=Count, y=reorder(Term, Count), fill=PValue))+
geom_bar(stat="identity", width=0.5)+
scale_fill_gradient(low="#FDB462", high="lightgrey")+
theme_bw()+
theme(text=element_text(family="Times"), axis.text.x=element_text(face="bold", size=14), axis.text.y=element_text(face="bold", size=14), axis.title.x=element_text(face="bold", size=15), axis.title.y=element_blank())






# Malaria-dsRNA-Transcriptome-Project
GO and KEGG analysis method

We primarily installed and used the “clusterProfiler”, "org.Pf.plasmo.db", “AnnotationHub”, "AnnotationDbi" packages in R to do the GO and KEGG analysis. “clusterProfiler” is a package that implements methods to analyze and visualize functional profilers (GO and KEGG) of gene and gene clusters [1].  Since Bioconductor has already provided OrgDb for about 20 species http://bioconductor.org/packages/release/BiocViews.html#___OrgDb, which made users conveniently build OrgDb via “AnnotationHub” [2]. We specifically used the “org.Pf.plasmo.db”, which includes genome wide annotation for Malaria. “AnnotationDbi” is an R package that provides an interface for connecting and querying various annotation databases [3]. 

After we input our data, we did the analysis on each of the three condition comparisons and named as “rnk_RS”, “rnk_RT”, and “rnk_TS”, which represents “Ring vs. Schizont”, “Ring vs. Troph”, and “Troph vs. Schizont”. We focused on the LogFC value of each condition comparisons and set the column names of the matrix as the GeneID. Then we removed the NA value if the GeneID and ranked it in descending order by fold change. Then we applied the gseGO() function, which represents Gene Set Enrichment Analysis of Gene Ontology, to get the GO analysis results under each of the condition comparisons. We set the p-value cut off threshold as 1 to receive a full results and to prevent some auto-screening of potentially important genes.    

Next, we did the KEGG analysis for each of the condition comparisons. clusterProfiler provides search_kegg_organisim() function to help searching supported organisms and the full list of KEGG supported organisms can be accessed via  http://www.genome.jp/kegg/catalog/org_list.html [4]. We used Plasmodium falciparum 3D7 in our analysis. We applied the gseKEGG() function, which represents Gene Set Enrichment Analysis of KEGG, to receive the KEGG analysis results under each of the condition comparison. Same as the GO analysis, we still set the p-vallue cutoff threshold as the maximum value 1 to in case we miss any of the potentially important results.

We can take the file “GO_RS_R” as an example to have a quick view of the GO analysis results. Each of the file includes 10 columns, which includes:
1.	GO term ID
2.	Description
3.	setSize: Number of genes in the gene set after filtering out those genes not in the expression dataset
4.	enrichmentScore: which reflects the degree to which a gene set is overrepresented at the top or bottom of a ranked list of genes. The ES is the maximum deviation from zero encountered in walking the list. A positive ES indicates gene set enrichment at the top of the ranked list; a negative ES indicates gene set enrichment at the bottom of the ranked list. [7]
5.	NES:  normalized enrichment score, that is, the enrichment score for the gene set after it has been normalized across analyzed gene sets. By normalizing the enrichment score, GSEA accounts for differences in gene set size and in correlations between gene sets and the expression dataset; therefore, the normalized enrichment scores (NES) can be used to compare analysis results across gene sets [7]. 
6.	pvalue
7.	p.adjust: not exactly equals to the qvalues [5]
8.	qvalues
9.	rank: The position in the ranked list at which the maximum enrichment score occurred. The more interesting gene sets achieve the maximum enrichment score near the top or bottom of the ranked list; that is, the rank at max is either very small or very large [7].
10. leading_edge: Displays the three statistics used to define the leading edge subset:
●     Tags. The percentage of gene hits before (for positive ES) or after (for negative ES) the peak in the running enrichment score. This gives an    indication of the percentage of genes contributing to the enrichment score.
●     List. The percentage of genes in the ranked gene list before (for positive ES) or after (for negative ES) the peak in the running enrichment score. This gives an indication of where in the list the enrichment score is attained.
●     Signal. The enrichment signal strength that combines the two previous statistics: where N is the number of genes in the list and Nh is the number of genes in the gene set. If the gene set is entirely within the first Nh positions in the list, then the signal strength is maximal or 100%. If the gene set is spread throughout the list, then the signal strength decreases towards 0%.
11.	core_enrichment: This is the subset of genes that contributes most to the enrichment result.

To take the file “KEGG_RS” as an example to go through the KEGG analysis results. It can be found that it looked similar as GO results, which has the same 11 columns as the previous file.  
Strictly speaking, the q-value and the FDR adjusted p-value are similar but not quite the same. However, the terms q-value and FDR adjusted p-value are often used generically by the Bioconductor community to refer to any quantity that controls or estimates any definition of the FDR.  In this general sense the terms are synonyms [6].

References:
[1]. https://www.bioconductor.org/packages//2.10/bioc/html/clusterProfiler.html
[2]. https://yulab-smu.github.io/clusterProfiler-book/chapter5.html
[3]. https://hbctraining.github.io/DGE_workshop_salmon/lessons/AnnotationDbi_lesson.html
[4]. https://yulab-smu.github.io/clusterProfiler-book/chapter6.html
[5]. https://www.biostars.org/p/128931/
[6]. https://stat.ethz.ch/pipermail/bioconductor/2012-December/049902.html
[7]. https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html



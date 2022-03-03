# Variant Analysis (DNA-seq): Tools Demonstration

Upon the request of Prof. Jacques-P. Tremblay, from the University of Laval, I prepared this demo to show my ability to analyze WGS data for discovering SNPs in genes. (May-Dec. 2016)
I demonstrate in this project the usage of 3 tools employed in the analysis of genetic variants:  
1. [variant_calling.R](variant_calling.R) is an R script that demonstrates using VariantTools R package which provides a mechanisim for calling variants.
The code is obtained from the tutorial written by Michael Lawrence  <http://bioconductor.org/packages/release/bioc/vignettes/VariantTools/inst/doc/tutorial.pdf>.

2. [variant_calling_samtools-mpileup.script](variant_calling_samtools-mpileup.script) is a bash script that demonstrates using samtools mpileup for obtaining a summary of the coverage of mapped reads on a reference sequence at a single base pair resolution, and BCFtools  which provides a mechanisim for calling variants. 
The code is obtained from Dave Tang’s blog <https://davetang.org/muse/2015/08/26/samtools-mpileup/>

3. [variant_annotation.R](variant_annotation.R) is an R script that demonstrates using VariantAnnotation R package which is used for annotating variants.  
The code for this file is obtained from the tutorial written by Valerie Obenchain <http://www.bioconductor.org/packages/release/bioc/vignettes/VariantAnnotation/inst/doc/VariantAnnotation.pdf>.

                           

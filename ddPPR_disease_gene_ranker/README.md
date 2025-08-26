This pipeline described an algorithm used to rank genes in disease based on experimental MS data. It expands an experimental interactome module based on a gene-associated-disease network to capture relevant genes that might me dysregulated or have causal relationships with a given phenotype in disease. The pipeline was originally used to compare interactome results of membrane WT- vs. F508del-CFTR to retrieve possible causal genes in known phenotypes of rescued F508del-CFTR, as well as possible dysregulated genes that might be associated with disease progress in CF, such as exacerbated inflammation. 

The pipeline consists of 1) generating a disease-associated network of protein-protein interactions that will be used for the downstream diffusion and ranking of genes and 2) the actual diffusion algorithm consisting of a comparative Personalized Page Rank that contrast two MS/MS interactome data for WT- and F508del-CFTR (Matos et al., 2019) as the seed set of each diffusion round for the same network, depicting which genes seem to be most differentially associated with one interactome relative to the other.

This project is under active development — for now it can only be used by cloning this repository.

Installation

Clone the repository:

git clone https://github.com/g-jp/projects/edit/main/ddPPR_disease_gene_ranker
cd ddPPR_disease_gene_ranker

Install localy:

pip install -e .

After installation, you can import functions directly:

from diffusion import PPR_diffusion, load_graph

References:

Matos, A. M., Pinto, F. R., Barros, P., Amaral, M. D., Pepperkok, R., & Matos, P. (2019). Inhibition of calpain 1 restores plasma membrane stability to pharmacologically rescued Phe508del-CFTR variant. Journal of Biological Chemistry, 294(36), 13396–13410. https://doi.org/10.1074/jbc.RA119.008738

READMEv1.0	REPO-NGS	https://github.com/ryfa5051/NGS.git

Next Generation Sequencing Project Files - which sadly solely amounts to lone bash scripts at the moment.

Including useful .sh files as well for things like FATSA manipulation for continued reference for my poor brain.

Current pipeline for Miseq reads generated courtesy of the NGS Core at JSCBB - Univerisity of Colorado at Boulder 

Changes:

*11-5-2018*
  - runtime now displays correctly in minutes
  - changed tirmmomatic back to default "false" for the <keepBothReads> reverse complement clone option
  - cleaned and annotated  
  
*done*

*11-14-2018*
  - added possible solution for misphasing with DNAsnp 
  - phasing is now done via Beagle5.0 and subsequent haplotype .vcf is manipulated via vcftools
  - this produces 2 FASTAs per individual containing allelic specific variants
  
 *done*
 

## Cufflinks class_code 
cuffClass=list(
	`Known`=c(
		`Complete match`="=",
		`Contained`="c"
	),
	`Novel`=c(
		`Potentially novel isoform`="j",
		`Within a reference intron`="i",
		`Exon on the opposite strand`="x",
		`Generic exonic overlap`="o",
		`Unkown, intergenic transcript`="u"
	),
	`Artefact`=c(
		`Multiple classifications`=".",
		`Possible polymerase run-on`="p",
		`Possible pre-mRNA fragment`="e",
		`Intron on the opposite strand`="s",
		`Repeat`="r"
	)
)
dt.cuffClass=data.table::data.table(
	rbind(
		  c(class="Known",class_code="=",desc="Complete match"),
		  c(class="Known",class_code="c",desc="Contained"),

		  c(class="Novel",class_code="j",desc="Potentially novel isoform"),
		  c(class="Novel",class_code="i",desc="Within a reference intron"),
		  c(class="Novel",class_code="x",desc="Exon on the opposite strand"),
		  c(class="Novel",class_code="o",desc="Generic exonic overlap"),
		  c(class="Novel",class_code="u",desc="Unkown, intergenic transcript"),

		  c(class="Artefact",class_code=".",desc="Multiple classifications"),
		  c(class="Artefact",class_code="p",desc="Possible polymerase run-on"),
		  c(class="Artefact",class_code="e",desc="Possible pre-mRNA fragment"),
		  c(class="Artefact",class_code="s",desc="Intron on the opposite strand"),
		  c(class="Artefact",class_code="r",desc="Repeat"))
)
dt.cuffClass$class=factor(dt.cuffClass$class,levels=c("Known","Novel","Artefact"))

#http://www.ensembl.org/common/Help/Glossary?db=core
#http://www.ensembl.org/Help/Faq?id=468
#http://www.ensembl.org/Help/View?id=151
bioType=list(
	`Protein coding`=c(
		'protein_coding', 'nonsense_mediated_decay', 'nontranslating_CDS', 'non_stop_decay', 'polymorphic_pseudogene', 'LRG_gene',
		'IG_C_gene', 'IG_D_gene', 'IG_gene', 'IG_J_gene', 'IG_LV_gene', 'IG_M_gene', 'IG_V_gene', 'IG_Z_gene', 
		'TR_C_gene', 'TR_D_gene', 'TR_J_gene', 'TR_V_gene'
	),
	`Pseudogene`=c(
		'pseudogene', 'processed_pseudogene', 'translated_processed_pseudogene', 'transcribed_processed_pseudogene', 'transcribed_unprocessed_pseudogene', 'unitary_pseudogene', 'transcribed_unitary_pseudogene', 
		'unprocessed_pseudogene', 'disrupted_domain', 'retained_intron',
		'IG_C_pseudogene', 'IG_D_pseudogene', 'IG_J_pseudogene', 'IG_V_pseudogene', 'IG_pseudogene',
		'TR_C_pseudogene', 'TR_D_pseudogene', 'TR_J_pseudogene', 'TR_V_pseudogene'
	),
	`Long noncoding`=c(
		'3prime_overlapping_ncrna', 'ambiguous_orf', 'antisense', 'antisense_RNA', 'lincRNA', 'ncrna_host', 
		'processed_transcript', 'sense_intronic', 'sense_overlapping', 'macro_lncRNA', 'vaultRNA'
	),
	`Short noncoding`=c(
		'miRNA', 'miRNA_pseudogene', 'piRNA', 'misc_RNA', 'misc_RNA_pseudogene', 'Mt_rRNA', 'Mt_tRNA', 'rRNA', 'scRNA', 'snlRNA', 
		'sRNA', 'snoRNA', 'scaRNA', 'snRNA', 'tRNA', 'tRNA_pseudogene', 'ribozyme'
	),
	`Others`=c(
		'TEC', 'Artifact'
	)
)

sampleFreq=c(`>=95%`=.95, `>=75%`=.75, `>=50%`=.5, `>=25%`=.25, `>=5%`=.05)

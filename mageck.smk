rule mageck_rra:
	input:	counts=os.path.join(output_dir,"counts_table/all_samples_assigned_count.csv")
	output:	genesummary=os.path.join(output_dir,"screen_analysis/{case}_vs_con.gene_summary.txt"),
		sgrnasummary=os.path.join(output_dir,"screen_analysis/{case}_vs_con.sgrna_summary.txt"),
		pdf=os.path.join(output_dir,"screen_analysis/{case}_vs_con.pdf")
	params:	prefix=os.path.join(output_dir,"screen_analysis/{case}_vs_con"),
		treatment=lambda wildcards: ",".join(CASE_IDS[wildcards.case]),
		control=lambda wildcards: ",".join(CONTROL_IDS),
		norm=config["norm_method"],
		controlsg=("" if not "neg_ctl" in config["library"] else "--control-sgrna "+config["library"]["neg_ctl"]),
		#cnv_correct=("" if not config["correct_cnv"] else "--cnv-norm "+config["cnv_norm"]+" --cell-line "+config["cnv_cell_line"]),
		additionalparameter=("" if not "additional_rra_parameter" in config else " "+config["additional_rra_parameter"])
	log:	"logs/mageck/{case}.rra.log"
	shell:
		"""
		module load statistical/R/3.5.3/gcc.8.3.1
		mageck test --norm-method {params.norm} \
			--output-prefix {params.prefix} \
			--count-table {input.counts} \
			{params.controlsg} \
			--treatment-id {params.treatment} \
			--control-id {params.control} {params.additionalparameter} 2> {log}
		"""

"""
Pipeline to process dual-CRISPR screen experiments
"""
import glob
import numpy
import os
import yaml
from Bio.Seq import Seq
from collections import defaultdict

configfile: "config.yaml"

# load cluster config file
CLUSTER = yaml.load(open(config['CLUSTER_YAML']), Loader=yaml.FullLoader)
wd = os.getcwd()
config["wd"] = wd

#Get relative position of Snakefile from wd
SNAKEFILE = workflow.snakefile
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)
config["snakefile"] = SNAKEFILE
config["snakefile_dir"] = SNAKEFILE_DIR

#Establish snakefile and environment dictionaries
rules_dir = os.path.abspath(os.path.join(SNAKEFILE_DIR, "rules"))               
scripts_dir = os.path.abspath(os.path.join(SNAKEFILE_DIR, "scripts"))                   
logs_dir = os.path.abspath(os.path.join(SNAKEFILE_DIR, "logs"))
output_dir = os.path.abspath(os.path.join(SNAKEFILE_DIR, "results"))

input_files = []

def postprocess_config(config):
	config["replicates"] = {}
	if "samples" in config:
		for sample in config["samples"]:
			replicates = config["samples"][sample]
			if not isinstance(replicates, list):
				replicates = [replicates]
			config["samples"][sample] = []
			for i, replicate in enumerate(replicates):
				name = "{}_{}".format(sample, i)
				config["replicates"][name] = {}
				config["replicates"][name]['r1'] = replicate
				config["samples"][sample].append(name)
	if "paired" in config:
		for sample in config["paired"]:
			paired_rep = config["paired"][sample]
			if not isinstance(paired_rep, list):
				paired_rep = [paired_rep]
			config["paired"][sample] = []
			for i, replicate in enumerate(paired_rep):
				name = "{}_{}".format(sample, i)
				config["replicates"][name]['r2'] = replicate
				config["paired"][sample].append(name)	
	# promoter U6: GATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCG
	# promoter H1: CGGGAAAGAGTGGTCTCATACAGAACTTAT, H1 reverse complement: ATAAGTTCTGTTTGAGACCACTCTTTCCCG
	# scaffold from 5'/3', GTTTTAGAGCTAGAAATAGCAAGTTGAGACGNNNNNNNNNNCGTCTCAACTTGCTATTTCTAGCTCTAAAAC
	config["pgrnas"]["P1fw"]=Seq(config["pgrnas"]["U6fw"])
	config["pgrnas"]["P1rv"]=config["pgrnas"]["P1fw"].reverse_complement()
	config["pgrnas"]["P2rv"]=Seq(config["pgrnas"]["H1rv"])
	config["pgrnas"]["P2fw"]=config["pgrnas"]["P2rv"].reverse_complement()

postprocess_config(config)
output_files = []
SAMPLE_IDS=list(config["replicates"].keys())
REPLICATE_IDS=list(config["paired"].keys())
CONTROL_IDS=config["experiments"]["rra"]["control"]
CASE_IDS={}
for id in REPLICATE_IDS:
        if id not in CONTROL_IDS and id not in config["experiments"]["rra"]["plasmid"]:
                if id.split('_rep')[0] not in CASE_IDS:
                        CASE_IDS[id.split('_rep')[0]]=[]
                        CASE_IDS[id.split('_rep')[0]].append(id)
                else:
                        CASE_IDS[id.split('_rep')[0]].append(id)
if "paired" in config:
        output_files.extend(expand(os.path.join(output_dir,"trimmed_reads/{sample}_R1.trim.fastq"), sample=SAMPLE_IDS))
        output_files.extend(expand(os.path.join(output_dir,"trimmed_reads/{sample}_R2.trim.fastq"), sample=SAMPLE_IDS))
#output_files.extend(expand(os.path.join(output_dir,"alignments/{sample}.sorted.bam"), sample=SAMPLE_IDS))
#output_files.extend(expand(os.path.join(output_dir,"alignments/{sample}.sorted.bam.bai"), sample=SAMPLE_IDS))
#output_files.extend(expand(os.path.join(output_dir,"alignments/{sample}.sorted.bam.flagstat"), sample=SAMPLE_IDS))
#output_files.extend(expand(os.path.join(output_dir,"alignments_filtered/{sample}.PEaligns_filtered.tsv"), sample=SAMPLE_IDS))
#output_files.extend(expand(os.path.join(output_dir,"alignments_filtered/{replicate}.reps_merged.tsv"), replicate=REPLICATE_IDS))
output_files.extend([os.path.join(output_dir,"counts_table/all_samples_assigned_count.csv"), os.path.join(output_dir,"counts_table/template_swapping_statistics.txt")])
#output_files.extend(expand(os.path.join(output_dir,"screen_analysis/{case}_vs_con.gene_summary.txt"), case=CASE_IDS))
#output_files.extend(expand(os.path.join(output_dir,"screen_analysis/{case}_vs_con.sgrna_summary.txt"), case=CASE_IDS))
#output_files.extend(expand(os.path.join(output_dir,"screen_analysis/{case}_vs_con.pdf"), case=CASE_IDS))

def get_dicts():
	dict_list = []
	for i in config['library']:
		if i.split('_')[-1] == 'dict':
			dict_list.append(config['library'][i])
	return dict_list


include: os.path.join(rules_dir, "mageck.smk")

rule all:
        input: output_files
        message: "Rule all"

if "samples" in config and "paired" in config:
	rule cutadapt_pe:
		input:  R1 = lambda wildcards: config['replicates'][wildcards.sample]['r1'],
			R2 = lambda wildcards: config['replicates'][wildcards.sample]['r2']
		output: R1p = os.path.join(output_dir,"trimmed_reads/{sample}_R1.trim.fastq"),
			R2p = os.path.join(output_dir,"trimmed_reads/{sample}_R2.trim.fastq")
		log:    os.path.join(logs_dir,"cutadapt_pe/{sample}.trim.log")
		threads: CLUSTER["cutadapt_pe"]["cpus-per-task"]
		shell:
			"""
			cutadapt -j {threads} -q 15,10 -e 0.1 -l {config[pgrnas][len]} --minimum-length=20 --pair-filter=any \
				-g {config[pgrnas][P1fw]} -g {config[pgrnas][P2rv]} -G {config[pgrnas][P2rv]} -G {config[pgrnas][P1fw]} \
				-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
				-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
				-A AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
				-a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
				-o {output.R1p} -p {output.R2p} {input.R1} {input.R2}
			"""
	rule bowtie2_map:
		input:	R1p = os.path.join(output_dir,"trimmed_reads/{sample}_R1.trim.fastq"),
			R2p = os.path.join(output_dir,"trimmed_reads/{sample}_R2.trim.fastq")
		output:	os.path.join(output_dir,"alignments/{sample}.sorted.bam")
		log:	os.path.join(logs_dir,"bowtie2_map/{sample}.align.log")
		threads: CLUSTER["bowtie2_map"]["cpus-per-task"]
		message: "align {input} to fake reference: {threads} threads"
		shell:
			"""
			bowtie2 --threads={threads} -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 --end-to-end \
				-x {config[library][bowtie2_index]} -1 {input.R1p} -2 {input.R2p} | samtools sort -O BAM -o {output} 2> {log}
			"""

	rule index_bam:
		input:  os.path.join(output_dir,"alignments/{sample}.sorted.bam")
		output: os.path.join(output_dir,"alignments/{sample}.sorted.bam.bai")
		log:    os.path.join(logs_dir,"index_bam/{sample}.index.log")
		threads: 1
		message:"index_bam {input}: {threads} threads"
		shell:
			"""
                	samtools index {input} > {output} 2> {log}
			"""

	rule flagstat_bam:
		input:  os.path.join(output_dir,"alignments/{sample}.sorted.bam")
		output: os.path.join(output_dir,"alignments/{sample}.sorted.bam.flagstat")
		log:    os.path.join(logs_dir,"flagstat_bam/{sample}.flagstat.log")
		threads: 1
		message:"flagstat_bam {input}: {threads} threads"
		shell:
			"""
			samtools flagstat {input} > {output} 2> {log}
			"""

	rule filter_bam:
		input:  os.path.join(output_dir,"alignments/{sample}.sorted.bam")
		output: os.path.join(output_dir,"alignments_filtered/{sample}.PEaligns_filtered.tsv")
		params:	tmp = os.path.join(output_dir,"alignments_filtered/{sample}.F2308f0x1.bam")
		log:    os.path.join(logs_dir,"filter_bam/{sample}.filter.log")
		threads: 1
		message:"filter_bam {input}: {threads} threads"
		shell:
			"""
			samtools view -F 2308 -f 1 -b {input} -o {params.tmp} 2> {log}
			join -1 1 -2 1 -a1 -a2 -e0 -o'0,1.3,1.5,1.6,1.7,1.9,1.10,2.3,2.5,2.6,2.7,2.9,2.10' <(samtools view -f 64 {params.tmp} | sort) <(samtools view -f 128 {params.tmp} | sort) > {output}
			rm {params.tmp}
			"""
	
rule merge_technical_replicates:
	input: lambda wildcards: [os.path.join(output_dir,"alignments_filtered/{0}.PEaligns_filtered.tsv".format(technical_id)) for technical_id in config["samples"][wildcards.replicate]]
	output: os.path.join(output_dir,"alignments_filtered/{replicate}.reps_merged.tsv")
	shell:
		"""
		cat {input} > {output}
		"""



rule generate_count_table:
	input:	counts = expand(os.path.join(output_dir,"alignments_filtered/{replicate}.reps_merged.tsv"), replicate=REPLICATE_IDS),
		dicts = get_dicts()
	output:	csv = os.path.join(output_dir,"counts_table/all_samples_assigned_count.csv"),
		stat = os.path.join(output_dir,"counts_table/template_swapping_statistics.txt")
	params:	qTHR = 2, QTHR = 40, mapQ = 1
	run:
		table, stat = {}, {}
		sample_ids=[]
		for countF in input.counts:
			iBase = os.path.basename(countF).split('.reps')[0]
			sample_ids.append(iBase)
			stat[iBase] = {}
			stat[iBase]['template'], stat[iBase]['switch'] = 0, 0
		for dictF in input.dicts:
			FIRSTLINE=True
			with open(str(dictF), 'r+') as d:
				for line in d:
					if FIRSTLINE:
						FIRSTLINE=False
						continue
					line=line.rstrip().split('\t')
					guide, region = str(line[0]), str(line[1])
					table[guide] = {}
					table[guide]['sgRNA'], table[guide]['gene'] = guide, region
					for s in sample_ids:
						table[guide][s]=0
		print("# Empty dictionary construction, DONE!")
		for countF in input.counts:
			with open(str(countF), 'r+') as f:
				s = os.path.basename(countF).split('.reps')[0]
				for line in f:
					line=line.rstrip().split()
					if line[1]==line[7]:
						stat[s]['template'] += 1
						if int(line[2])>=int(params.qTHR) and int(line[8])>=int(params.qTHR): 
							table[str(line[1])][s] += 1
					else:
						if line[1]!='0' and line[7]!='0':
							if int(line[2])>=int(params.mapQ) and int(line[8])>=int(params.mapQ):
								stat[s]['switch'] += 1
			print("# Count hits assign for sample: {}, DONE!".format(s))
		colnames=['sgRNA','gene']+sample_ids
		print(','.join(colnames), file=open(str(output.csv), 'a+'))
		for g in table:
			print(','.join(str(table[g][c]) for c in colnames), file=open(str(output.csv), 'a+'))
		print('\t'.join(['Samples', '#template', '#switch']), file=open(str(output.stat), 'a+'))
		for s in stat:
			print('\t'.join([str(s), str(stat[s]['template']), str(stat[s]['switch'])]), file=open(str(output.stat), 'a+'))

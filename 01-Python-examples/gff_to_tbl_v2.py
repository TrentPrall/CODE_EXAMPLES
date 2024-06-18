from BCBio import GFF
import argparse
import csv

"""
This script converts a gff file of genomic annotations output by Geneious to NCBI tbl format for
easy submission of large annotation files into genbank. This is to accompany submission of reference
fasta contigs.

The script takes two arguments:
1. the input gff file
2. a product table corresponding to the gene annotations. This is necessary for tbl format.
"""
# Set up argument parsing
parser = argparse.ArgumentParser(description="Convert gff to NCBI tbl format for genbank submission")
parser.add_argument("in_file", type=str, help="Path to the input gff file.")
parser.add_argument("product_table", type=str, help="Path to the tsv product table.")

args = parser.parse_args()

in_file = args.in_file
product_table = args.product_table

# append protein products
product_dict = {}

with open(product_table, "r") as product_object:
	product_reader = csv.reader(product_object, delimiter = "\t")
	for line in product_reader:
		product_dict[line[0]] = line[1]
product_object.close()

def make_feature_per_sequence(fasta_header):
	# {Name:[type, start, Stop, Strand, product, [CDS intervals] }
	# establish ordered gene, pseudogene, and ncRNA dictionaries
	gff_dict = {}

	# parse input gff file, collect metadata, and organize it into python dict
	in_handle = open(in_file)
	parse = GFF.parse(in_handle)
	for rec in parse:
		for feature in rec.features:
			if feature.type == "promoter" or feature.type == "regulatory" or feature.type == "motif" or feature.type == "exon":
				pass
			else:
				# print(str(feature.location).split(']')[0].strip('[').split(':')[0])
				type = feature.type
				#print(type)
				Name = feature.qualifiers['Name'][0]
				Start_parse = int(str(feature.location).split(']')[0].strip('[').split(':')[0])
				#print(Start_parse)
				Stop = int(str(feature.location).split(']')[0].strip('[').split(':')[1])
				Strand = str(feature.location).split('](')[-1].strip(')')
				#print(Strand)
				product = product_dict[Name]
				
				if "KIR" in Name:
					Name = Name.split('_')[0]
				else:
					Name = Name
				
				
				#print(rec.id)
				# iterate only through ids that match the input sequence name
				if rec.id == fasta_header:
					Start = Start_parse + 1
					if type == "gene":
						gff_dict[Name] = ["gene",Start, Stop, Strand, product, []]
			
					elif type == "pseudogene":
						gff_dict[Name] = ["pseudo",Start, Stop, Strand, product, []]
			
					elif type == "pseudo_special":
						gff_dict[Name] = ["pseudo",Start, Stop, Strand, product, []]
			
					elif type == "ncRNA_gene":
						gff_dict[Name] = ["ncRNA",Start, Stop, Strand, product, []]
				else:
					pass
		
					
	in_handle.close()
	return gff_dict
	
	
def append_CDS_intervals(fasta_header,gff_dict):
	in_handle = open(in_file)
	parse = GFF.parse(in_handle)
	# append CDS intervals
	for rec in parse:
		for feature in rec.features:
			type = feature.type
			name = feature.qualifiers['Name'][0]
			#print(name)
			interval1 = ''
			
			if type == "CDS":
				if "KIR" in name:
					name = name.split('_')[0]
				else:
					name = name
				
				if rec.id == fasta_header:	
					print(name)
					interval1 = str(feature.location).strip('[)').split('](')[0]
					gff_dict[name][5].append(interval1)
	in_handle.close()

		
			
	

### BEGIN MAIN EXECUTION
in_handle = open(in_file)
parse = GFF.parse(in_handle)

headers = []
for rec in parse:
	headers.append(rec.id)

in_handle.close()
print(headers)

# open output file
output_object = open("output.tbl", 'w')
output_writer = csv.writer(output_object, delimiter = '\t')

for sequence in headers:
	print(sequence)
	fasta_gff_dict = make_feature_per_sequence(sequence)
	print(fasta_gff_dict)
	append_CDS_intervals(sequence, fasta_gff_dict)
	

	# Sort Dictionary based on start postion
	sortedbyvalue = {k: v for k, v in sorted(fasta_gff_dict.items(), key=lambda v: v[1][2])}
	# print(sortedbyvalue)
	
	# WRITE TBL OUTPUT
	sequence_name = [">Feature " + sequence]
	output_writer.writerow(sequence_name)
	
	for key in sortedbyvalue:
		name = key
		type = sortedbyvalue[key][0]
		start1 = sortedbyvalue[key][1]
		stop1 = sortedbyvalue[key][2]
		strand = sortedbyvalue[key][3]
		product = sortedbyvalue[key][4]
		exons = sortedbyvalue[key][5]
		
		
		# SORT CDS INTERVALS BY POSITION
		sorted_intervals = {}
		for interval in exons:
			start = int(interval.split(':')[0])
			stop = int(interval.split(':')[1])
			sorted_intervals[start] = stop
		sortedintervalsbystart = {k: v for k, v in sorted(sorted_intervals.items())}
		
		# RULES FOR GENES WITH CDS ENTERIES
		if type == "gene":
			if strand == "+":
				outline = [str(start1), str(stop1), 'gene','',''] # write gene interval
				output_writer.writerow(outline)
				
				outline = ['','','',"gene", name]
				output_writer.writerow(outline) # write gene name identifier
				
				# write out CDS intervals
				iter = 0
				for interval in sortedintervalsbystart:
					if  iter == 0: # write out first interval with CDS denoted
						outline = [interval + 1, sortedintervalsbystart[interval], "CDS",'','']
						output_writer.writerow(outline)
						iter = iter + 1
					else:
						outline = outline = [interval + 1, sortedintervalsbystart[interval], '','','']
						output_writer.writerow(outline)
						
				outline = ['','','',"gene", name]	
				output_writer.writerow(outline) # write gene name identifier
				
				outline = ['','','',"product", product]	
				output_writer.writerow(outline) # write gene name identifier
				
				
			elif strand == "-":
				outline = [str(stop1), str(start1), 'gene','',''] # write gene interval
				output_writer.writerow(outline)
	
				outline = ['','','',"gene", name]
				output_writer.writerow(outline) # write gene name identifier
	
				# write out CDS intervals
				iter = 0
				for interval in exons:
					if  iter == 0: # write out first interval with CDS denoted
						outline = [interval.split(':')[1], int(interval.split(':')[0]) + 1, "CDS",'','']
						output_writer.writerow(outline)
						iter = iter + 1
					else:
						outline = [interval.split(':')[1], int(interval.split(':')[0]) + 1,'','','']
						output_writer.writerow(outline)
		
				outline = ['','','',"gene", name]	
				output_writer.writerow(outline) # write gene name identifier
	
				outline = ['','','',"product", product]	
				output_writer.writerow(outline) # write gene name identifier
	
		# RULES FOR PSEUDOGENES
		elif type == "pseudo":
			if strand == "+":
				outline = [str(start1), str(stop1), 'gene','',''] # write gene interval
				output_writer.writerow(outline)
			elif strand == "-":
				outline = [str(stop1), str(start1), 'gene','',''] # write gene interval
				output_writer.writerow(outline)
	
			outline = ['','','',"gene", name]	
			output_writer.writerow(outline) # write gene name identifier
	
			outline = ['','','',"gene_desc", product]	# substitute gene_desc for product only for pseudogene entries 
			output_writer.writerow(outline) # write gene name identifier
	
			outline = ['','','', "pseudo", '']
			output_writer.writerow(outline) # write gene name identifier
		
		
		
		# RULES FOR NCRNA
		elif type == "snoRNA" or type == "miRNA" or type =="lncRNA":
			if strand == "+":
				outline = [str(start1), str(stop1), 'gene','',''] # write gene interval
				output_writer.writerow(outline)
	
				outline = ['','','',"gene", name]
				output_writer.writerow(outline) # write gene name identifier
	
				# write out CDS intervals
				iter = 0
				for interval in sortedintervalsbystart:
					if  iter == 0: # write out first interval with CDS denoted
						outline = [interval + 1, sortedintervalsbystart[interval], "ncRNA",'','']
						output_writer.writerow(outline)
						iter = iter + 1
					else:
						outline = outline = [interval + 1, sortedintervalsbystart[interval], '','','']
						output_writer.writerow(outline)
		
				outline = ['','','',"ncRNA_class", type]	
				output_writer.writerow(outline) # write gene name identifier		
				
				outline = ['','','',"gene", name]	
				output_writer.writerow(outline) # write gene name identifier
	
				outline = ['','','',"product", product]	
				output_writer.writerow(outline) # write gene name identifier
	
	
			elif strand == "-":
				outline = [str(stop1), str(start1), 'gene','',''] # write gene interval
				output_writer.writerow(outline)
	
				outline = ['','','',"gene", name]
				output_writer.writerow(outline) # write gene name identifier
	
				# write out CDS intervals
				iter = 0
				for interval in exons:
					if  iter == 0: # write out first interval with CDS denoted
						outline = [interval.split(':')[1], int(interval.split(':')[0]) + 1, "ncRNA",'','']
						output_writer.writerow(outline)
						iter = iter + 1
					else:
						outline = [interval.split(':')[1], int(interval.split(':')[0]) + 1,'','','']
						output_writer.writerow(outline)
	
	
				outline = ['','','',"ncRNA_class", type]	
				output_writer.writerow(outline) # write gene name identifier		
	
				outline = ['','','',"gene", name]	
				output_writer.writerow(outline) # write gene name identifier
	
				outline = ['','','',"product", product]	
				output_writer.writerow(outline) # write gene name identifier
output_object.close()
	



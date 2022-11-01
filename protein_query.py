#!/usr/bin/python3
"""
CREATED BY USER: B195515-2021
PROGRAM: Python for protein sequence analysis
FOR BPSM ICA2
"""

""" Setup python3 environment, functions, shortcuts """	
import os, sys, subprocess, re
import numpy as np
import pandas as pd
from glob import glob

def shell(cmd):
	""" Executes a shell command """
	return subprocess.call(cmd, shell = True)

def shellout(cmd):
	""" Saves output of shell command """
	return subprocess.check_output(cmd, shell = True).decode("utf-8").rstrip('\n')

def convert_spaces(var, character):
	""" Removes spaces in string, joined by 'character' """
	if ' ' in var:
		var = (character).join(var.split(' '))
		return var
	else:
		return var

def show_output(distype, filename, extension):
	""" Takes user input to display output of a section 
	distype: name of program to display output eg. 'display' or 'firefox' (str)
	filename: base name of the file (str)
	extension: extension of the file without the first dot (str)
	"""
	while True:
		input_display = input('Display output? Yes/No\n\t')
		if input_display in ans_yes:
			cmd_display = f'{distype} {filename}.{extension}'
			shell(cmd_display)
			print(f'\nOutput saved as {filename}.{extension}\n')
			break
		elif input_display in ans_no:
			print(f'\nOutput saved as {filename}.{extension}\n')
			break
		else:
			print(f'Invalid input: {input_display}. Please try again.\n')
			continue

def exit_message():	
	files = glob('*')
	files.sort(key=os.path.getmtime)
	print(f'{lines()}\nEND.\nThank you for using this program.\nFollowing are the output files generated.\n{lines()}')
	print('Directory: ', os.getcwd())
	print('\n'.join(files))
	print('\nType \'exit()\' to exit python3 or press enter to continue in interactive shell.')
	return

def lines():
	return '-'*70

ans_yes = {'Yes','yes','YES','yup','yep','OK','ok','yeah','yea','1','y','Y'}
ans_no = {'No','NO','no','nope','nah','nuh uh','exit','0','n','N'}
ans_symbol = r'[,\\({\)[\}\].?;:"\'`!@#$%&*~/^]'

""" Make new directory for this analysis """
print('Welcome to the Python program for protein sequence analysis by B195515-2021.\n')
while True: 
	dirname = convert_spaces(input("Enter name for output directory (eg. my_protein_analysis1):\n\t"),'_')
	if dirname in os.listdir(): 
		# must create new directory
		print(f'Directory exists \'{dirname}\'. Input a different directory name.\n')
		continue
	elif ( re.search(ans_symbol, dirname) ) or dirname == '': 
		# cannot be empty/invalid chars
		print(f'Invalid input {dirname}. Please try again.\n')
		continue
	else:
		# creates directory and change workspace
		os.mkdir(dirname)
		os.chdir(dirname)
		print(f'{lines()}\nThank you. You are currently here:\n{lines()}\n\t{os.getcwd()}\n\t{dirname}\n\nAll outputs for this run will be saved in this folder for your reference.')
		break

""" Initial inputs """
while True: 
	taxgrp = input("Enter taxon group (eg. birds):\n\t")
	protfam = input("Enter protein family (eg. glucose-6-phosphatase):\n\t")
	if ( re.search(ans_symbol, taxgrp) ) or ( re.search(ans_symbol, protfam) ) or taxgrp == '' or protfam == '': # cannot be empty/invalid chars
		print(f'Invalid input:\t{taxgrp} and\t{protfam}. Please try again.\n')
		continue
	else: # create a base name for all output files
		print(f'{lines()}\nThank you. Input has been received. You have chosen:\nProtein family: {protfam}\nTaxon/taxID: {taxgrp}\n{lines()}')
		basename = convert_spaces(taxgrp,'-')+'_'+convert_spaces(protfam,'-')
		break

""" Determine partial/non-partial query """
print('Now using EDIRECT to retrieve sequences from the NCBI database for your query\n')
while True: 
	input_partial = input("Input 'partial' to include partial sequences or 'NOT partial' for only full sequences:\n\t")
	if input_partial == 'partial': # get all sequences
		print(f'{lines()}\nThank you. All sequences will be searched for.\n{lines()}')
		input_partial = ''
		break
	elif input_partial == 'NOT partial': # get non partial sequences
		print(f'{lines()}\nThank you. Only non partial sequences will be searched for.\n{lines()}')
		break
	else: # catches all other inputs
		print(f'Invalid input: {input_partial}. Please try again.\n')
		continue

""" Generate query fasta file """
print(f'Now generating fasta file and accession ID file for query {taxgrp} and {protfam}.\n\nPlease be patient :)\n')
shell(fr'esearch -db protein -query "{taxgrp}[Organism:exp] AND {protfam}*[Protein] NOT isoform NOT predicted {input_partial}"| efetch -format fasta > {basename}.prot.fa')

try: 
	# Check if a file is created AND check that it's not empty
	if os.stat(f'{basename}.prot.fa').st_size != 0:
		shell(fr'esearch -db protein -query "{taxgrp}[Organism:exp] AND {protfam}*[Protein] NOT isoform NOT predicted {input_partial}"| efetch -format acc > {basename}.prot.acc')
		print('Successful! Fasta file and accession ID file containing all matching sequences from your query has been created in this directory\n')
		cmd_seqfound = int(shellout('grep -c ">" *.prot.fa'))
		print(f'{lines()}\nThere are {cmd_seqfound} matches to your query\n{lines()}')
	else:
		sys.exit('Query failure. No match found.\nRe-run program with a valid taxon group or protein family.')
except:
	sys.exit('Query failure. No match found.\nRe-run program with a valid taxon group or protein family.')
	

""" Split multiseq fastafile into list of single sequences with headers """
singleseq = open(f'{basename}.prot.fa').read().rstrip('\n').split('>')[1:]

""" Create query dictionary """
singleseq_dict = {}
for i, entry in enumerate(singleseq):
	if entry == []:
		print('No sequence found')
		break
	else:
		entry = '>'+entry.replace('\n','')
	if entry.startswith(">sp") or entry.startswith(">tr"):
		print(f'Skipping sequence {i} retrieved from Swiss-Prot and TrEMBL database.\n')
		continue
	if entry.startswith(">"):
		# Dict Key
		seq_id = re.search(r'^>\S{0,40}',entry).group().replace('>','')
		if seq_id not in singleseq_dict:
			# Create Key in Dict
			singleseq_dict[seq_id] = []
		# Dict Values
		species = re.search(r'(?<=\[).+?(?=\])', entry).group()
		pattern = fr'{seq_id}\s(.*?)\s\['
		prot_name = re.search(pattern, entry).group(1)
		pattern2 = fr'\[{species}\]'
		sequence = re.search('(.+)'+pattern2+'(.+)', entry).group(2)
		# Add Values to the dictionary
		singleseq_dict[seq_id].append(prot_name)	
		singleseq_dict[seq_id].append(species)
		singleseq_dict[seq_id].append(sequence)

"""Create query dataframe"""
print('Creating query dataframe...\n')
df = pd.DataFrame.from_dict(singleseq_dict, orient='index',
	columns =['Protein_name', 'Species', 'AA_sequence'])
# set acc id as column, not index
df = df.rename_axis('Accession_ID').reset_index()

""" Save query dataframe to csv file """
try:
	while True:
		input_savedf = input('Save original dataframe to a tab-separated csv file? Yes(recommended)/No\n\t')
		if input_savedf in ans_yes:
			df.to_csv(f'{basename}.dataframe.csv', sep='\t', na_rep='NA')
			print(f'\n{basename}.dataframe.csv file created.\n')
			show_output('firefox', basename, 'dataframe.csv')
			break
		elif input_savedf in ans_no:
			print('Dataframe not saved to file.')
			break
		else:
			print(f'\nInvalid input: {input_savedf}. Please try again.\n')
			continue
except:
	print('ERROR in saving dataframe.\n')


"""Show matches to query"""
species_names = set(df['Species'].unique())
print('Following are the summary matches to your query...')
print(f'{lines()}\nSummary\n{lines()}\n',df.describe())
print(f'{lines()}\nTop 15 protein matches\n{lines()}\n',df['Protein_name'].value_counts()[0:15])
print(f'{lines()}\nFollowing are the {len(species_names)} species found in your query\n{lines()}\n', df['Species'].value_counts().to_string())

def check_species_input(var):
	""" Checks if excluded species are in the species set"""
	count = 0
	for i in var:
		if i not in species_names:
			print(f'\nInvalid input {i}\n')
			count += 1
		else:
			count += 0
	if count > 0:
		return 'invalid'
	else:
		return 'good'

""" Let user choose species """
while True:
	species_exclude = input('\nEnter species name(s) to EXCLUDE from analysis, separated by comma(s).\nElse, enter \'NONE\' to include all species:\n\t')
	if species_exclude == '': # reject empty input
		print('\nNo input. Please try again.\n')
		continue
	elif species_exclude.upper() == 'NONE': # FINAL: use all species
		print(f'\nIncluding all {len(species_names)} species in the analysis.\n')
		species_final = species_names
		break
	else: # Split by commas, remove trailing spaces, capitalize first letter and get unique values
		species_exclude = set(i.strip().capitalize() for i in species_exclude.split(','))
	if check_species_input(species_exclude) == 'invalid': # re-enter input
		continue
	else: # FINAL: subset species
		print(f'\nExcluding species {species_exclude} in the analysis.\n')
		species_final = list(species_names.difference(species_exclude))
		""" Reconstruct df and fasta sequences from chosen species """
		df2 = df[df['Species'].isin(species_final)]
		df2.to_csv(f'{basename}.dataframe2.csv', sep='\t', na_rep='NA')
		print(f'\n{basename}.dataframe2.csv file created.\n')
		show_output('firefox', basename, 'dataframe2.csv')
		singleseq2  = []
		for seq in singleseq:
			species = re.search(r'(?<=\[).+?(?=\])', seq).group()
			if species in species_final:
				singleseq2.append(seq)
			else:
				pass
		singleseq = singleseq2
		break

""" Let user choose number of sequences """
length = len(singleseq)
if length < 2:
	print('Only one sequence found. Nothing to align.\n')
	exit_message()
if 2 <= length <= 1000:
	input_q = input(f'{length} matches found. Trim number of queries? Yes/No\n\t')
	if input_q in ans_yes:
		try:
			input_qlength= int(input('Input integer value (Max 1000):\n\t'))
			singleseq = singleseq[:input_qlength]
		except:
			print(f'All {length} queries will be processed')
	else: print(f'All {length} queries will be processed')
if 1000 < length <= 4000: # Limit up to 4000
	input_q = input(f'{length} matches found. Trim number of queries? Yes/No\n\t')
	if input_q in ans_no: # do not trim sequences
			print(f'All {length} queries will be processed\n')
	if input_q in ans_yes: # Limit to 1000 sequences
		try:
			input_qlength= int(input('Input integer value (Max 4000):\n\t'))
			singleseq = singleseq[:input_qlength]
		except:
			print('Processing 1000 queries instead.\n')
			singleseq = singleseq[:1000]
	else:
		print('Invalid input.\n')
		print('Processing 1000 queries instead.\n')
		singleseq = singleseq[:1000]
if length > 4000:
	input_q = input(f'{length} matches found. Trim required.\nInput integer value (Max 4000, Recommended 1000):\n\t')
	try:
		if int(input_q) <= 4000:
			print(f'Query of {input_q} matches will be processed.\n')
			singleseq = singleseq[:int(input_q)]
	except:
		print('Invalid input.\n')
		print('Processing 1000 queries instead.\n')
		singleseq = singleseq[:1000]

""" Create new fasta/accession file if trimmed """
for seq in singleseq:	
	with open(f'{basename}.protfinal.fa','a') as finalseq:
		finalseq.write(f'>{seq}')
	with open(f'{basename}.protfinal.acc', 'a') as finalseqacc:
		finalseqacc.write(f'{seq.split()[0]}\n')
print(f'{lines()}\nTwo new files created:\n\t{basename}.protfinal.fa\n\t{basename}.protfinal.acc\n{lines()}')

""" Multiple sequence alignment with Clustalo """
try:
	print("\nGenerating alignment of fasta protein sequences...\n")
	shell(fr'clustalo -i {basename}.protfinal.fa -v -t protein --guidetree-out={basename}.guidetree --percent-id --threads=10 -o {basename}.clust.fa --force')
	if f'{basename}.clust.fa' in os.listdir():
		print('\nSuccessful! Query sequence alignment file created.\n')
except:
	print('ERROR creating sequence alignment.')
	sys.exit('Exiting the program. Please try again.')

""" Wildcard analysis #1: EMBOSS Infoalign"""
print(f'{lines()}\nInfoalign displays basic information about a multiple sequence alignment\n{lines()}')
while True:
	input_infoalign = input('\nRun Infoalign analysis? Yes/No\n\t')
	if input_infoalign in ans_yes:
		print('\nRunning infoalign on alignment and output as html file.\n')
		shell(fr'infoalign {basename}.clust.fa {basename}.infoalign -html -auto')
		print('Infoalign run successful.\n')
		show_output('firefox', basename, 'infoalign')
		break
	elif input_infoalign in ans_no:
		print('Skippping infoalign analysis.\n')
		break
	else:
		print(f'Invalid input: {input_infoalign}. Please try again.\n')
		continue

""" Plot conservation of amino acids with plotcon """	
print(f'{lines()}\nNow plotting for amino acid sequence conservation for your query.\n{lines()}\n')
while True:
	plottitle = str(input("Input plotcon image file name. eg. birds_glucose-6-phosphatase\n\t"))
	if ( re.search(ans_symbol, plottitle) ):	# non-symbols
		print(f'Invalid input: {plottitle}. Please try again.\n')
		continue
	else:
		plottitle = convert_spaces(plottitle,'-')	# remove spaces
		break		
while True:
	try:
		winsize = int(input("Input an integer for window size. (Default[4], larger= smoother plot, Max=10)\n\t"))
		if winsize < 2 or winsize > 10:	# integer val between 2 and 10
			print(f'Invalid input value of {winsize}. Please try again.')
			continue
		else:
			break
	except:
		print('Invalid non-integer input. Please try again.')	# integer input only
		continue

print(f'{lines()}\nInput received\n{lines()}\nProgram name:\tEMBOSS Plotcon\nFile name:\t{basename}.clust.fa\nWindow size:\t{winsize}\nOutput file:\t{plottitle}\nDefault matrix used:\tBLOSUM62\n')
shell(f"plotcon {basename}.clust.fa -winsize {winsize} -goutfile {plottitle} -graph cps")
show_output('display', plottitle, 'ps')


""" Identify protein motifs using PROSITE(prosextract) and patmatmotifs """
print(f'{lines()}\nNow searching for protein motifs using prosextract and patmatmotifs\n{lines()}\n')

""" Import prosite database """
cmd_prolines = r'cp /localdisk/software/EMBOSS-6.6.0/share/EMBOSS/data/PROSITE/prosite.lines .'
prosurl = "https://ftp.expasy.org/databases/prosite/"
shell(fr'wget -qO prosite.dat {prosurl}')
shell(fr'wget -qO prosite.doc {prosurl}')
if 'prosite.dat' in os.listdir() and 'prosite.doc' in os.listdir():
	print(f'FTP file prosite.dat and prosite.doc successfully downloaded from {prosurl} to current folder.\n')
	#shell(r'prosextract . -auto -stdout True -warning True -version True') # unhash if prosite.lines need to be made
	shell(cmd_prolines)
else: 
	print(f'FTP from {prosurl} is unsuccessful.\nUsing pre-loaded database from MSC5...\n')
	shell(cmd_prolines)

""" Get input to prune sequences before using patmatmotifs """
print('Pruning in patmatmotifs ignores all common post-translational modifications:\nmyristyl\nasn_glycosylation\ncamp_phospho_site\npkc_phospho_site\nck2_phospho_site\ntyr_phospho_site.\n')
while True:
	input_prune = input("Prune basic motifs from analysis? 'Yes'(recommended)/No:\n\t")
	if input_prune in ans_yes:
		print('\nBasic motifs will not be included in the output file.\n')
		input_prune = 'Yes'
		break
	elif input_prune in ans_no:
		print('\nBasic motifs will be included in the output file.\n')
		input_prune = 'No'
		break
	else:
		print(f'Invalid input: {input_prune}. Please try again.\n')
		continue
print('Thank you.\n\nProceeding to process motif recognition...\n\nPlease be patient :)\n')

""" Run patmatmotifs """
for i,seq in enumerate(singleseq):
	with open('singleseq.temp','w') as out1:	# Create a temp input file
		out1.write(f'>{seq}')
	shell(fr'patmatmotifs singleseq.temp motif.temp -prune {input_prune} -auto Yes')	# Execute patmatmotif
	with open('motif.temp') as out3:	
		singlemotif = out3.read().rstrip('\n')
		nomotifs = singlemotif.find('HitCount: 0') # Classify sequence with/without motifs; irrelevant if not pruned
	if nomotifs == -1:
		with open(f'{basename}_with.motifs','a') as out4:
			out4.write(singlemotif)
	else:
		with open(f'{basename}_without.motifs','a') as out5:
			out5.write(singlemotif)

if i+1 == len(singleseq):	# check for end of sequence processing
	print(f'Patmatmotif ran on all {i+1} sequences.\n\nTwo output files are created:\n\t{basename}_with.motifs\n\t{basename}_without.motifs')
else:
	print(f'Patmatmotif ran on {i+1} sequences. Error processing some sequences.\n')



""" Output the number of unique motifs found in all sequences """
cmd_motifcount = fr'grep "Motif = " {basename}_with.motifs| sort| uniq -c'
print(f'{lines()}\nUnique amino acid motifs found\n{lines()}\n',shellout(cmd_motifcount))
while True:
	input_motifcount = input('\nCreate .motifcount file? Yes(recommended)/No\n\t')
	if input_motifcount in ans_yes:
		with open(f'{basename}.motifcount','w') as out4:
			out4.write('The number of motifs NOT from post-translational modification from your query is as follows:\n\nTaxon group:\t'+taxgrp+'\nProtein family:\t'+protfam+'\n\nCount\tMotif\n'+shellout(cmd_motifcount))
			print('\nA .motifcount file created.\n')
			show_output('firefox', basename, 'motifcount')
			break
	elif input_motifcount in ans_no:
		print('No .motifcount file created\n')
		break
	else:
		print('Invalid input detected:\t', input_motifcount)
		continue


""" Wildcard analysis #2: EMBOSS Pepnet """
print(f'{lines()}\nPepnet draws an amino acid helical net for a protein sequence.\n{lines()}\n')
while True:
	input_runpepnet = input('Run pepnet for your query and save image? Yes/No\n\t')
	if input_runpepnet in ans_yes:
		if f'{basename}.clust.fa' in os.listdir():
			print('\nGenerating helical net representation for the alignment.\n')
			shell(fr'pepnet {basename}.clust.fa -goutfile {basename}.pepnet -graph cps -auto')
			show_output('display', basename, 'pepnet.ps')
			break
		else:
			print('\nAlignment file not found.\nSkipping pepnet analysis.\n')
			break
	elif input_runpepnet in ans_no:
		print('\nSkipping pepnet analysis.\n')
		break
	else:
		print('\nInvalid input detected:\t', input_runpepnet)
		continue

""" Wildcard analysis #3: EMBOSS pepwheel """
print(f'{lines()}\nPepwheel plots the properties of alpha helices in a protein,\nillustrating the orientation of each amino acid according to its polarity.\nThis structure may reveal patterns indicating a specific motif or protein fold.\n{lines()}\n')
while True:
	input_runpepwheel = input('Run pepwheel for your query and save image? Yes/No:\n\t')
	if input_runpepwheel in ans_yes:
		if f'{basename}.clust.fa' in os.listdir():
			print('\nGenerating helical wheel diagram for the alignment.\n')
			shell(fr'pepwheel {basename}.clust.fa -gtitle Helical_wheel -gsubtitle {basename} -goutfile {basename}.pepwheel -graph cps -auto')
			show_output('display', basename, 'pepwheel.ps')
			break
		else:
			print('\nAlignment file not found.\nSkipping pepwheel analysis.\n')
			break
	elif input_runpepwheel in ans_no:
		print('\nSkipping pepwheel analysis.\n')
		break
	else:
		print('\nInvalid input detected:\t', input_runpepwheel)
		continue

""" Wildcard analysis #4: EMBOSS pepwindowall """
print(f'{lines()}\nPepwindowall produces a Kyte-Doolittle hydropathy plot\nfrom an aligned set of protein sequences\n{lines()}\n')
while True:
	input_runpepwindow = input('Run pepwindowall for your query and save image? Yes/No:\n\t')
	if input_runpepwindow in ans_yes:
		if f'{basename}.clust.fa' in os.listdir():
			print('\nGenerating hydropathy plot for the alignment.\n')
			shell(fr'pepwindowall {basename}.clust.fa -gtitle Hydropathy_plot -gsubtitle {basename} -goutfile {basename}.pepwindow -graph cps -auto')
			show_output('display', basename, 'pepwindow.ps')
			break
		else:
			print('\nAlignment file not found.\nSkipping pepwindow analysis.\n')
			break
	elif input_runpepwindow in ans_no:
		print('\nSkipping pepwindow analysis.\n')
		break
	else:
		print('\nInvalid input detected:\t', input_runpepwindow)
		continue

""" Wildcard analysis #5: EMBOSS antigenic """
print(f'{lines()}\nAntigenic predicts the potential antigenic site(s)\nof the homologous sequence of a protein family\n{lines()}\n')
while True:
	input_runant = input('Run antigenic for your query and save image? Yes/No:\n\t')
	if input_runant in ans_yes:
		if f'{basename}.clust.fa' in os.listdir():
			print('\nGenerating an antigen \'motif\' report for the alignment.\n')
			shell(fr'antigenic {basename}.clust.fa {basename}.antigenic -sprotein Yes -auto')
			show_output('firefox', basename, 'antigenic')
			break
		else:
			print('\nAlignment file not found.\nSkipping antigenic analysis.\n')
			break
	elif input_runant in ans_no:
		print('\nSkipping antigenic analysis.\n')
		break
	else:
		print('\nInvalid input detected:\t', input_runant)
		continue

""" Cleanup """
shell('rm -fr *.temp prosite*')
""" Exit """
exit_message()

"""END OF SCRIPT"""
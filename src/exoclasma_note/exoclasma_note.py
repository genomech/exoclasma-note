__scriptname__ = 'exoclasma-note'
__version__ = 'v0.9.0'

import json
import logging
import multiprocessing
import os
import pandas #
import subprocess
import tempfile

# -----=====| LOGGING |=====-----

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.DEBUG)

# ------======| MISC |======------

def MultipleTags(Tag, List, Quoted = True):
	Result = list()
	for Item in List:
		Result.append(Tag)
		Result.append(ArmorDoubleQuotes(Item) if Quoted else Item)
	return ' '.join(Result)

def ArmorDoubleQuotes(String): return f'"{String}"'

def ArmorSingleQuotes(String): return f"'{String}'"

# ------======| SUBPROCESS |======------

def BashSubprocess(SuccessMessage, Command, AllowedCodes = list()):
	logging.debug(f'Shell command: {Command}')
	Shell = subprocess.Popen(Command, shell = True, executable = 'bash', stdout = subprocess.PIPE, stderr = subprocess.PIPE)
	_, Stderr = Shell.communicate()
	if (Shell.returncode != 0) and (Shell.returncode not in AllowedCodes):
		logging.error(f'Shell command returned non-zero exit code {Shell.returncode}: {Command}\n{Stderr.decode("utf-8")}')
		exit(1)
	if (Shell.returncode in AllowedCodes):
		logging.warning(f'Shell command returned ALLOWED non-zero exit code {Shell.returncode}')
	logging.info(SuccessMessage)

# -----=====| ANNOVAR |=====-----

def ANNOVAR(InputVCF, OutputVCF, AnnovarFolder, GenomeAssembly, Threads = multiprocessing.cpu_count()):
	with tempfile.TemporaryDirectory() as TempDir:
		TableAnnovarPath = os.path.join(AnnovarFolder, 'table_annovar.pl')
		AnnovarDBPath = os.path.join(AnnovarFolder, 'humandb')
		TempVCF = os.path.join(TempDir, 'temp.vcf')
		AnnotatedVCF = f'{TempVCF}.{GenomeAssembly}_multianno.vcf'
		# Compose commands
		TempVcfCommand = ['zcat', ArmorDoubleQuotes(InputVCF), '>', ArmorDoubleQuotes(TempVCF)]
		AnnotationCommand = ['perl', ArmorDoubleQuotes(TableAnnovarPath), ArmorDoubleQuotes(TempVCF), 
				ArmorDoubleQuotes(AnnovarDBPath), '--buildver', GenomeAssembly, '--protocol', 'refGene,ensGene,knownGene',
				'--operation', 'g,g,g', '--remove', '--vcfinput', '--thread', str(Threads)]
		BgzipCommand = ['bgzip', '-c', ArmorDoubleQuotes(AnnotatedVCF), '>', ArmorDoubleQuotes(OutputVCF)]
		TabixCommand = ['tabix', '-p', 'vcf', ArmorDoubleQuotes(OutputVCF)]
		# Processing
		BashSubprocess('VCF file copied', ' '.join(TempVcfCommand))
		BashSubprocess('ANNOVAR annotation finished', ' '.join(AnnotationCommand), AllowedCodes = [25])
		BashSubprocess('VCF bgzipped', ' '.join(BgzipCommand))
		BashSubprocess('VCF indexed', ' '.join(TabixCommand))

def AnnovarStage(Unit):
	ANNOVAR(
		InputVCF = os.path.join(Unit['OutputDir'], Unit['Output']['VCF']),
		OutputVCF = os.path.join(Unit['OutputDir'], Unit['Output']['AnnovarVCF']),
		AnnovarFolder = Unit['AnnovarFolder'],
		GenomeAssembly = Unit['Reference']['GenomeInfo']['annovar_alias'],
		Threads = Unit['Config']['Threads']
		)

def AnnotationPipeline(UnitFile, AnnovarFolder, AnnovarGenomeAlias):
	logging.info(f'{__scriptname__} Annotate {__version__}')
	ConfigPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.json')
	UnitPath = os.path.realpath(UnitFile)
	Config = json.load(open(ConfigPath, 'rt'))
	Unit = json.load(open(UnitPath, 'rt'))
	logging.info(f'Unit loaded: "{UnitPath}"')
	Unit['Output']['AnnovarVCF'] = f'_temp.{Unit["ID"]}.annovar.vcf.gz'
	StageAlias = 'Annovar'
	if StageAlias not in Unit['Stage']:
		Unit['AnnovarFolder'] = os.path.realpath(AnnovarFolder)
		Unit['Reference']['GenomeInfo']['annovar_alias'] = str(AnnovarGenomeAlias)
		AnnovarStage(Unit)
		Unit['Stage'].append(StageAlias)
		json.dump(Unit, open(UnitFile, 'wt'), indent = 4, ensure_ascii = False)
	StageAlias = 'Tabix'
	if StageAlias not in Unit['Stage']:
		
	logging.info('Job finished')

AnnotationPipeline('/mydocs/MyDocs/Cloud/Core/exoclasma-note/tests/testdata/EmilTestExoC/unit.json', '/home/thousandtowers/annovar/annovar', 'hg19')

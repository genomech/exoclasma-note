__scriptname__ = 'exoclasma-note'
__version__ = 'v0.9.0'

import gzip
import json
import logging
import multiprocessing
import os
import re
import subprocess
import tempfile

SamplesFunc = {
	'DP': lambda x: int(x),
	'AD': lambda x: [int(k) for k in x.split(',')],
	'PL': lambda x: [int(k) for k in x.split(',')],
	'GT': lambda x: [int(k) for k in re.split('[\|/]', x)],
	'PGT': lambda x: [int(k) for k in re.split('[\|/]', x)],
	'PID': lambda x: str(x),
	'GQ': lambda x: int(x),
	'PS': lambda x: int(x),
	}

VcfInfoFunc = {
	'AC': lambda x: [int(k) for k in x.split(',')],
	'AF': lambda x: [float(k) for k in x.split(',')],
	'AN': lambda x: int(x),
	'BaseQRankSum': lambda x: float(x),
	'DP': lambda x: int(x),
	'ExcessHet': lambda x: float(x),
	'FS': lambda x: float(x),
	'MLEAC': lambda x: [int(k) for k in x.split(',')],
	'MLEAF': lambda x: [float(k) for k in x.split(',')],
	'MQ': lambda x: float(x),
	'MQRankSum': lambda x: float(x),
	'QD': lambda x: float(x),
	'ReadPosRankSum': lambda x: float(x),
	'SOR': lambda x: float(x)
	}

def GeneDetailFunc(Line):
	if Line == '.': return None
	if 'dist\\x3d' not in Line: return Line.split('\\x3b')
	Result = [k.split('\\x3d') for k in Line.split('\\x3b')]
	Result = [int(k[1]) if k[1] != 'NONE' else None for k in Result]
	return Result

def AAChangeFunc(Line):
	if (Line == '.') or (Line == 'UNKNOWN'): return None
	List = Line.split(',')
	Result = list()
	for Item in List:
		Record = dict()
		Splitted = Item.split(':')
		assert len(Splitted) == 5
		Record['gene'] = str(Splitted[0])
		Record['accession'] = str(Splitted[1])
		assert Splitted[2][:4] == 'exon'
		Record['exon'] = int(Splitted[2][4:])
		assert Splitted[3][:2] == 'c.'
		Record['transcript'] = str(Splitted[3])
		assert Splitted[4][:2] == 'p.'
		Record['protein'] = str(Splitted[4])
		Result.append(Record)
	return Result

AnnovarFunc = {
	'Func': lambda x: x.split('\\x3b'),
	'Gene': lambda x: [(k if k != 'NONE' else None) for k in x.split('\\x3b')],
	'GeneDetail': GeneDetailFunc,
	'ExonicFunc': lambda x: None if x == '.' else str(x),
	'AAChange': AAChangeFunc
	}

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

def GetStream(File):
	Stream = gzip.open(File, 'rt')
	while True:
		Line = next(Stream)
		if Line[:2] != '##': break
	return Stream

def ExtractSamples(Vcf):
	Stream = gzip.open(Vcf, 'rt')
	while 1:
		Row = next(Stream)
		if Row[:6] == "#CHROM": break
	return {index: item for index, item in enumerate(Row[:-1].split('\t')[9:])}

def ParseVcfRow(Row, Samples):
	try:
		Qual = int(Row[5])
	except ValueError:
		try:
			Qual = float(Row[5])
		except ValueError:
			Qual = None
	try:
		Format = Row[8]
	except IndexError:
		Format = None
	else:
		if Format == '.': Format = None
	Result = {
		"CHROM": Row[0],
		"POS": int(Row[1]),
		"ID": None if (Row[2] == '.') else Row[2],
		"REF": Row[3],
		"ALT": Row[4].split(','),
		"QUAL": Qual,
		"FILTER": None if (Row[6] == '.') else ([] if (Row[6] == 'PASS') else Row[6].split(';')),
		"INFO": Row[7]
	}
	INFO = [i.split('=') for i in Result["INFO"].split(';')]
	Result["INFO"] = {}
	Result["ANNOTATIONS"] = []
	Stage = 0
	for item in INFO:
		if item[0] == 'ANNOVAR_DATE':
			Stage += 1
			Result["ANNOTATIONS"].append({'ANNOVAR': { 'refGene': {}, 'ensGene': {}, 'knownGene': {}}})
			continue
		if item[0] == 'ALLELE_END': continue
		if Stage == 0: Result["INFO"][item[0]] = VcfInfoFunc[item[0]](item[1])
		else:
			TheWay = item[0].split('.')
			Result["ANNOTATIONS"][-1]['ANNOVAR'][TheWay[1]][TheWay[0]] = AnnovarFunc[TheWay[0]](item[1])
	Format = Row[8].split(':')
	Result["SAMPLES"] = { Samples[index]: { Format[i]: SamplesFunc[Format[i]](value) for i, value in enumerate(item.split(':'))} for index, item in enumerate(Row[9:]) }
	return Result

def ConversionStage(Unit):
	Vcf = os.path.join(Unit['OutputDir'], Unit['Output']['AnnovarVCF'])
	OutFile = os.path.join(Unit['OutputDir'], Unit['Output']['VariantsJSON'])
	Samples = ExtractSamples(Vcf)
	Stream = GetStream(Vcf)
	Output = gzip.open(OutFile, 'wt')
	while True:
		try:
			Line = next(Stream)[:-1].split('\t')
		except StopIteration:
			break
		Record = ParseVcfRow(Line, Samples)
		Output.write(json.dumps(Record, separators = (',', ':'), ensure_ascii = False) + '\n')
	Stream.close()
	Output.close()

def AnnotationPipeline(UnitFile, AnnovarFolder, AnnovarGenomeAlias):
	logging.info(f'{__scriptname__} Annotate {__version__}')
	ConfigPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.json')
	UnitPath = os.path.realpath(UnitFile)
	Config = json.load(open(ConfigPath, 'rt'))
	Unit = json.load(open(UnitPath, 'rt'))
	logging.info(f'Unit loaded: "{UnitPath}"')
	Unit['Output']['AnnovarVCF'] = f'_temp.{Unit["ID"]}.annovar.vcf.gz'
	Unit['Output']['VariantsJSON'] = f'_temp.{Unit["ID"]}.variants.json.gz'
	StageAlias = 'Annovar'
	if StageAlias not in Unit['Stage']:
		Unit['AnnovarFolder'] = os.path.realpath(AnnovarFolder)
		Unit['Reference']['GenomeInfo']['annovar_alias'] = str(AnnovarGenomeAlias)
		AnnovarStage(Unit)
		Unit['Stage'].append(StageAlias)
		json.dump(Unit, open(UnitFile, 'wt'), indent = 4, ensure_ascii = False)
	StageAlias = 'Conversion'
	if StageAlias not in Unit['Stage']:
		ConversionStage(Unit)
	logging.info('Job finished')

AnnotationPipeline('/mydocs/MyDocs/Cloud/Core/exoclasma-note/tests/testdata/EmilTestExoC/unit.json', '/home/thousandtowers/annovar/annovar', 'hg19')

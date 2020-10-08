__author__ = 'Aitor Blanco (aitor.blancomiguez@unitn.it)'
__version__ = '0.01'
__date__ = '12 Aug 2020'

import os, re
import subprocess as sb
from utils import error

UNIREF_VERSION = '201906'

"""
Executes a command

:param cmd: the command to execute
"""
def execute(cmd):
    inp_f = None
    out_f = sb.DEVNULL

    if cmd['stdin']:
        inp_f = open(cmd['stdin'], 'r')
    if cmd['stdout']:
        out_f = open(cmd['stdout'], 'w')
    print(cmd['command_line'])

    sb.run(cmd['command_line'], stdin=inp_f, stdout=out_f)

    if cmd['stdin']:
        inp_f.close()
    if cmd['stdout']:
        out_f.close()


"""
Execute mmseq2.

:param config: the mmseq2 params
"""
def execute_mmseq2(config):
    params = {
        "program_name" : "mmseqs",
        "params" : config,
        "command_line" : "#program_name# #params#"
    }
    execute(compose_command(params))


"""
Execute prokka.

"""
def execute_prokka(input_file, outdir, prefix, nprocs=4):
    config = '--centre X --compliant --force'
    params = {
        "program_name" : "prokka",
        "input" : "",
        "output" : "--prefix",
        "output_path" : "--outdir",
        "threads" : "--cpus",
        "params" : config,
        "command_line" : "#program_name# #output_path# #output# #params# #threads# #input#"
    }
    execute(compose_command(params=params, input_file=input_file, output_path=outdir,
        output_file=prefix, nproc=nprocs))


"""
Execute uniref_annotator.

"""
def execute_uniref_annotator(input_file, output, tmpdir, DATABASES_PATH, DIAMOND_PATH, nprocs=4):
    config = '--seqtype prot --diamond {} --uniref90db {}/UniRef90_{}.dmnd --uniref50db {}/UniRef50_{}.dmnd --diamond-options "--threads {}"'.format(DIAMOND_PATH, DATABASES_PATH, UNIREF_VERSION, DATABASES_PATH, UNIREF_VERSION, nprocs)
    params = {
        "program_name" : os.path.join(os.path.dirname(os.path.realpath(__file__)), 'uniref_annotator', 'uniref_annotator.py'),
        "input" : "",
        "output" : "--out",
        "output_path" : "--temp",
        "params" : config,
        "command_line" : "#program_name# #output# #output_path# #params# #input#"
    }
    execute(compose_command(params=params, input_file=input_file, output_path=tmpdir,
        output_file=output, nproc=nprocs))


"""
Compress to tar.bz2

:param input: the input to compress
:param output_dir: the output directory
:returns: the compressed file
"""
def compress_tar_bz2(input, output_dir):
    params = {
        "program_name" : "tar",
        "params" : "-cjf",
        "database" : "",
        "command_line" : "#program_name# #params# #output# #database# #input#"
    }
    execute(compose_command(params, input_file=os.path.basename(input),
        output_file=os.path.join(output_dir, os.path.basename(input) + ".tar.bz2"),
        database="--directory=" + output_dir))
    return os.path.join(output_dir, os.path.basename(input) + ".tar.bz2")


"""
Execute bowtie2-build.

:param input: the input fasta file(s)
:param output: the output directory
"""
def execute_bowtie2_build(input, output, nprocs, BOWTIE2_PATH):
    params = {
        "program_name" : "bowtie2-build",
        "params" : "--threads " + str(nprocs),
        "command_line" : "#program_name# #params# #input# #output#"
    }
    execute(compose_command(params, input_file=input, output_file=output))


"""
Execute bowtie2-inspect.

:param input: the input fasta file(s)
:param output: the output directory
"""
def execute_bowtie2_inspect(input, BOWTIE2_PATH):
    params = {
        "program_name" : "bowtie2-inspect",
        "params" : "-n",
        "command_line" : "#program_name# #params# #input#"
    }
    execute(compose_command(params, input_file=input))


"""
Compress to BZ2

:param input: the file to compress
:returns: the compressed file
"""
def compress_bz2(input):
    params = {
        "program_name" : "bzip2",
        "command_line" : "#program_name# #input#"
    }
    execute(compose_command(params, input_file=input))


"""
Compose a command for further executions

:param params: the params of the command
:param check: [default=False] check if program is available
:param sub_mod: [optional] the model
:param input_file: [optional] the input file
:param database: [optional] the database
:param output_path: [optional] the output path
:param output_file: [optional] the output file
:param nproc: [default=1] the number of procs to use
"""
def compose_command(params, check=False, sub_mod=None, input_file=None, database=None, output_path=None, output_file=None, nproc=1):
    program_name = None
    stdin = None
    stdout = None
    environment = os.environ.copy()
    r_output_path = None
    r_output_file = None
    command_line = params['command_line']

    if 'program_name' in list(params):
        command_line = command_line.replace('#program_name#', params['program_name'])
        program_name = params['program_name']
    else:
        error('Error: something wrong... '+program_name+' not found!', exit=True)

    if check:
        command_line = program_name

        if 'version' in list(params):
            command_line = '{} {}'.format(program_name, params['version'])
    else:
        if 'params' in list(params):
            command_line = command_line.replace('#params#', params['params'])

        if 'threads' in list(params):
            command_line = command_line.replace('#threads#', '{} {}'.format(params['threads'], nproc))

        if output_path:
            r_output_path = output_path

            if 'output_path' in list(params):
                command_line = command_line.replace('#output_path#', '{} {}'.format(params['output_path'], output_path))
            else:
                output_file = os.path.join(output_path, output_file)

        if sub_mod:
            mod = sub_mod

            if 'model' in list(params):
                mod = '{} {}'.format(params['model'], sub_mod)

            command_line = command_line.replace('#model#', mod)

        if input_file:
            inp = input_file

            if 'input' in list(params):
                inp = '{} {}'.format(params['input'], input_file)

            if '<' in command_line:
                command_line = command_line.replace('<', '')
                command_line = command_line.replace('#input#', '')
                stdin = inp
            else:
                command_line = command_line.replace('#input#', inp)

        if database and ('database' in list(params)):
            command_line = command_line.replace('#database#', '{} {}'.format(params['database'], database))

        if output_file:
            out = output_file
            r_output_file = output_file

            if 'output' in list(params):
                out = '{} {}'.format(params['output'], output_file)

            if '>' in command_line:
                command_line = command_line.replace('>', '')
                command_line = command_line.replace('#output#', '')
                stdout = out
            else:
                command_line = command_line.replace('#output#', out)

        if 'environment' in list(params):
            new_environment = dict([(var.strip(), val.strip())
                                    for var, val in [a.strip().split('=') for a in params['environment'].split(',')]])
            environment.update(new_environment)

    # find string sourrunded with " and make them as one string
    quotes = [j for j, e in enumerate(command_line) if e == '"']

    for s, e in zip(quotes[0::2], quotes[1::2]):
        command_line = command_line.replace(command_line[s + 1:e], command_line[s + 1:e].replace(' ', '#'))

    return {'command_line': [str(a).replace('#', ' ') for a in re.sub(' +', ' ', command_line.replace('"', '')).split(' ') if a],
            'stdin': stdin, 'stdout': stdout, 'env': environment, 'output_path': r_output_path, 'output_file': r_output_file}

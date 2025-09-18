""" File handling configuration and setup for the snakemake pipeline

The important parameters from the config file are read in, and
the different wells and cycles are detected from rawinput/ or input.
Helper functions and wildcard constraints are also declared here
"""

import os
import sys
import glob
import time


# Read in directories from config file
rawinput_dir = config.get('rawinput_dir', 'rawinput/')
input_dir = config.get('input_dir', 'input/')
stitching_dir = config.get('stitching_dir', 'stitching/')
segmentation_dir = config.get('segmentation_dir', 'segmentation/')
sequencing_dir = config.get('sequencing_dir', 'sequencing/')
phenotyping_dir = config.get('phenotyping_dir', 'phenotyping/')
output_dir = config.get('output_dir', 'output/')
qc_dir = config.get('qc_dir', output_dir + 'qc/')

# The prefix used in rawinput/ for phenotype cycles
phenotype_date = config.get('phenotype_date', 'phenotype')
#phenotype_cycle = config.get('phenotype_cycle', 'PT')

phenotype_scale = config['phenotype_scale']
bases_scale = config['bases_scale']

if os.path.exists(rawinput_dir):
    print ('rawinput')
    dates = sorted(os.listdir(rawinput_dir))
    dates_pt = dates.copy()

    phenotype_dates = [date for date in dates if date[:len(phenotype_date)] == phenotype_date]
    if 'phenotype_cycles' not in config:
        phenotype_cycles = ['PT', 'PT1', 'PT2', 'PT3', 'PT4'][:len(phenotype_dates)]
    else:
        phenotype_cycles = config['phenotype_cycles']


    if 'wells' not in config:
        #wells = sorted([path.replace('Well', '').partition('_')[0] for path in os.listdir(rawinput_dir + '/' + dates[0])])
        wells = sorted(list(set([path.partition('Well')[2].partition('_')[0] for path in glob.glob(rawinput_dir + '/*/*.nd2')])))
    else:
        print (config['wells'])
        wells = config['wells']

    for date in phenotype_dates:
        if date in dates_pt:
            dates.remove(date)

    if 'cycles' not in config:
        cycles = ['{:02}'.format(i) for i in range(len(dates))]
    else:
        cycles = config['cycles']

else:
    if 'wells' not in config:
        wells = [dirname.replace('well', '') for dirname in sorted(os.listdir(input_dir)) if dirname != 'auxdata']

        for i in range(len(wells)):
            well = wells[i]
            well = well.split('_section')[0]
            well = well.split('_cyclenoise')[0].split('_noise')[0]
            well = well.split('_subset')[0]
            wells[i] = well

        wells = list(set(wells))
    else:
        wells = config['wells']

    if 'cycles' not in config:
        cycles_pt = [dirname[5:] for dirname in sorted(os.listdir(input_dir + '/well' + wells[0]))]
        cycles = [cycle for cycle in cycles_pt if cycle[0] != 'P']
        phenotype_cycles = [cycle for cycle in cycles_pt if cycle[0] == 'P']

        dates_pt = []#['date' + cycle for cycle in cycles_pt]
        dates = []#['date' + cycle for cycle in cycles]
        phenotype_dates = []
    else:
        cycles = config['cycles']
        phenotype_cycles = config['phenotype_cycles']

cycles_pt = cycles + phenotype_cycles
#cycles_pt = sorted(cycles_pt)

cellpose_cyto_index = config.get('cellpose_cyto_index', 1)
cellpose_diameter = config.get('cellpose_diameter', 50)
cellpose_cycle = config.get('cellpose_cycle', cycles[-1] if len(cycles) else None)

wildcard_constraints:
    well = '(well)?(' + '|'.join(wells) + ')(_subset\d+)?(_(cycle|)noise\d+)?(_section\d+)?',
    well_stitching = '(well)?(' + '|'.join(wells) + ')(_subset\d+)?(_(cycle|)noise\d+)?',
    well_nonoise = '(well)?(' + '|'.join(wells) + ')(_subset\d+)?',
    well_nosubset = '(well)?(' + '|'.join(wells) + ')',
    well_base = '(well)?(' + '|'.join(wells) + ')',

    tile = '\d\d\d\d',
    cycle = '|'.join(cycles_pt),

    path = '([^/]*/)*[^/.]*',
    path_nogrid = '((?!_grid\d)[^.])*',

    segmentation_type = 'cells|nuclei|cellsbases|nucleibases',

# Regex for use with wildcard constraints
phenotyping_channel_regex = '(' + '|'.join('{}|{}'.format(i, re.escape(name))
            for i, name in enumerate(config['phenotyping_channels'])) + ')'
sequencing_channel_regex = '(' + '|'.join('{}|{}'.format(i, name)
            for i, name in enumerate(config['sequencing_channels'])) + ')'
# only matches channels in both sequencing and phenotyping
any_channel_regex = ('(' + '|'.join(map(str, range(min(len(config['sequencing_channels']), len(config['sequencing_channels'])))))
    + '|' + '|'.join(map(re.escape, set(config['sequencing_channels']) & set(config['phenotyping_channels']))) + ')')

assert all(let in config['sequencing_channels'] for let in 'GTAC')

# slice that will extract all sequencing channels from sequencing images
sequencing_channels_slice = [config['sequencing_channels'].index(let) for let in 'GTAC']
if max(sequencing_channels_slice) - min(sequencing_channels_slice) == 3:
    sequencing_channels_slice = slice(min(sequencing_channels_slice), max(sequencing_channels_slice) + 1)

# order of sequencing channels
sequencing_channels_order = [chan for chan in config['sequencing_channels'] if chan in 'GTAC']


def debug(*args, **kwargs):
    print (time.asctime() + ':', *args, **kwargs, file=sys.stderr)

def progress(*args, **kwargs):
    import starcall.utils
    return starcall.utils.simple_progress(*args, **kwargs)

def print_mem(name, mem_limit):
    import psutil
    import starcall.utils
    #{key: starcall.utils.human_readable(val) for key,val in psutil.Process().memory_info()._asdict().items()}
    mem = psutil.Process().memory_info().rss / 1000000
    debug ("Rule", name, "Limit", mem_limit, starcall.utils.human_readable(mem_limit * 1000000),
                    "Using", mem, starcall.utils.human_readable(mem * 1000000))

def parse_param(name, default_value):
    def func(wildcards):
        val = getattr(wildcards, name)
        if val == '':
            return default_value
        val = val[len(name)+1:]
        try: return int(val)
        except: pass
        try: return float(val)
        except: pass
        return val
    return func

def param_constraint(name, pattern):
    return '(_' + name + '(' + pattern + '))?'

def params_regex(*params):
    return '(' + ''.join('(_{}[^_]+)?'.format(name) for name in params) + ')'

def channel_index(channel, kind=None, cycle=None):
    if kind is None and cycle is not None:
        kind = 'phenotyping' if cycle in phenotype_cycles else 'sequencing'
    if type(channel) == str:
        return config[kind+'_channels'].index(channel)
    return channel

def coredump():
    if os.fork() == 0:
        os.abort()

def print_info():
    print ()
    print ('FISSEQ data pipeline summary')
    print ('  rawinput directory:', rawinput_dir)
    print ('  input directory:', input_dir)
    print ('  stitching directory:', stitching_dir)
    print ('  sequencing directory:', sequencing_dir)
    print ('  phenotyping directory:', phenotyping_dir)
    print ('  output directory:', output_dir)
    print ()
    print ('Input found:')
    print ('  Cycles:')
    for cycle,date in zip(cycles_pt, dates_pt):
        print ('    cycle' + cycle + ':', rawinput_dir + date + '/')
    print ()
    print ('  Wells:', ', '.join(wells))
    print ()

rule print_info:
    run:
        print_info()

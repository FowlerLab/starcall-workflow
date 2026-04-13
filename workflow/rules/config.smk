""" File handling configuration and setup for the snakemake pipeline

The important parameters from the config file are read in, and
the different wells and cycles are detected from rawinput/ or input.
Helper functions and wildcard constraints are also declared here
"""

import os
import sys
import glob
import time
import matplotlib
matplotlib.use('agg')


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

##### Parsing alternative config options ####

if 'wells' in config and type(config['wells']) == int:
    config['wells'] = ['well{}'.format(i) for i in range(1, config['wells'] + 1)]
    print (config['wells'])

if 'cycles' in config:
    if type(config['cycles']) == int:
        config['cycles'] = ['{:02}'.format(i) for i in range(config['cycles'])]
    if type(config['cycles'][0]) == int:
        cycles = ['{:02}'.format(i) for i in config['cycles']]

if 'phenotype_cycles' in config:
    if type(config['phenotype_cycles']) == int:
        config['phenotype_cycles'] = (['PT'] + ['P{}'.format(i) for i in range(1, config['phenotype_cycles'])])[:config['phenotype_cycles']]


##### Finding all input files #####

if 'inputfiles' not in config:
    possible_files_input = []
    for root, dirs, files in os.walk(input_dir):
        possible_files_input.extend(root + '/' + filename for filename in files)

    possible_files_raw = []
    for root, dirs, files in os.walk(rawinput_dir, followlinks=True):
        possible_files_raw.extend(root + '/' + filename for filename in files)

    #possible_files = sorted(possible_files_input) + sorted(possible_files_raw)
    possible_files = sorted(possible_files_raw)

    detect_wells = 'wells' not in config
    inputfiles = {}

    for path in possible_files:
        if not (path.endswith('.tif') or path.endswith('.tiff') or path.endswith('.nd2')):
            continue

        if detect_wells:
            if path.count('well') or path.count('Well'):
                well = path.replace('Well', 'well')
                well = well[well.index('well'):].split('_')[0]
                inputfiles.setdefault(well, []).append(path)
        else:
            matching_wells = []
            for well in config['wells']:
                if any(option in path for option in [well, well[:4].replace('well', 'Well') + well[4:], 'well' + well, 'Well' + well]):
                    matching_wells.append(well)
            matching_wells = sorted(matching_wells, key=lambda path: len(path))
            if len(matching_wells) != 0:
                inputfiles.setdefault(matching_wells[-1], []).append(path)

    if detect_wells:
        config['wells'] = list(inputfiles.keys())

    detected_cycles = set()
    detected_pt_cycles = set()

    for well in inputfiles.keys():
        cyclepaths = {}
        index = 0
        pt_index = 0
        for path in inputfiles[well]:
            if phenotype_date in path:
                if 'phenotype_cycles' in config:
                    if pt_index >= len(config['phenotype_cycles']): continue
                    cycle = config['phenotype_cycles'][pt_index]
                else:
                    cycle = 'PT' if pt_index == 0 else 'P{}'.format(pt_index)
                pt_index += 1
                detected_pt_cycles.add(cycle)
            else:
                if 'cycles' in config:
                    if index >= len(config['cycles']): continue
                    cycle = config['cycles'][index]
                else:
                    cycle = '{:02}'.format(index)
                index += 1
                detected_cycles.add(cycle)

            cyclepaths[cycle] = path

        inputfiles[well] = cyclepaths

    if 'cycles' not in config:
        config['cycles'] = sorted(detected_cycles)
    if 'phenotype_cycles' not in config:
        config['phenotype_cycles'] = sorted(detected_pt_cycles)

    #for well, files in inputfiles.items():
        #print (well)
        #print ('\n'.join('\t{}: {}'.format(*pair) for pair in files.items()))

    #print (config['cycles'])
    #print (config['phenotype_cycles'])

    config['inputfiles'] = inputfiles


"""
if os.path.exists(rawinput_dir):
    dates = sorted(os.listdir(rawinput_dir))
    dates_pt = dates.copy()

    phenotype_dates = [date for date in dates if date[:len(phenotype_date)] == phenotype_date]
    if 'phenotype_cycles' not in config:
        phenotype_cycles = ['PT', 'P1', 'P2', 'P3', 'P4'][:len(phenotype_dates)]
    else:
        phenotype_cycles = config['phenotype_cycles']
        if type(phenotype_cycles) == int:
            phenotype_cycles = ['PT', 'P1', 'P2', 'P3', 'P4'][:phenotype_cycles]


    if 'wells' not in config:
        #wells = sorted([path.replace('Well', '').partition('_')[0] for path in os.listdir(rawinput_dir + '/' + dates[0])])
        wells = [os.path.basename(path) for path in glob.glob(rawinput_dir + '/*/*.nd2')]
        wells = [well.split('.nd2')[0].split('_Chan')[0] for well in wells]
        wells = sorted(list(set(wells)))
        #wells = sorted(list(set([os.path.basename(path).partition('_')[0] for path in glob.glob(rawinput_dir + '/*/*.nd2')])))
        wells = [well[0].replace('W', 'w') + well[1:] for well in wells]
    else:
        wells = config['wells']
        if type(wells) == int:
            wells = ['well{}'.format(well) for well in range(1, wells + 1)]

    for date in phenotype_dates:
        if date in dates_pt:
            dates.remove(date)

    if 'cycles' not in config:
        cycles = ['{:02}'.format(i) for i in range(len(dates))]
    else:
        cycles = config['cycles']
        if type(cycles) == int:
            cycles = ['{:02}'.format(i) for i in range(cycles)]
        if len(cycles) and type(cycles[0]) == int:
            cycles = ['{:02}'.format(i) for i in cycles]

else:
    if 'wells' not in config:
        wells = [dirname for dirname in sorted(os.listdir(input_dir)) if dirname != 'auxdata']

        for i in range(len(wells)):
            well = wells[i]
            well = well.split('_section')[0]
            well = well.split('_cyclenoise')[0].split('_noise')[0]
            well = well.split('_subset')[0]
            wells[i] = well

        wells = list(set(wells))
    else:
        wells = config['wells']
        if type(wells) == int:
            wells = ['well{}'.format(well) for well in range(1, wells + 1)]

    if 'cycles' not in config:
        first_well = glob.glob(input_dir + wells[0] + '*/')[0]
        cycles_pt = [dirname[5:] for dirname in sorted(os.listdir(first_well))]
        cycles = [cycle for cycle in cycles_pt if cycle[0] != 'P']
        phenotype_cycles = [cycle for cycle in cycles_pt if cycle[0] == 'P']

        dates_pt = []#['date' + cycle for cycle in cycles_pt]
        dates = []#['date' + cycle for cycle in cycles]
        phenotype_dates = []
    else:
        cycles = config['cycles']
        phenotype_cycles = config['phenotype_cycles']
        if type(cycles) == int:
            cycles = ['{:02}'.format(i) for i in range(cycles)]
        if type(cycles[0]) == int:
            cycles = ['{:02}'.format(i) for i in cycles]
        if type(phenotype_cycles) == int:
            phenotype_cycles = ['PT', 'P1', 'P2', 'P3', 'P4'][:phenotype_cycles]
"""

cycles = config['cycles']
phenotype_cycles = config['phenotype_cycles']
wells = config['wells']

cycles_pt = cycles + phenotype_cycles
#cycles_pt = sorted(cycles_pt)

cellpose_cyto_index = config.get('cellpose_cyto_index', 1)
cellpose_diameter = config.get('cellpose_diameter', 50)
cellpose_cycle = config.get('cellpose_cycle', cycles[-1] if len(cycles) else None)

wildcard_constraints:
    well = '(' + '|'.join(wells) + ')(_subset\d+)?(_split\d+)?(_(cycle|)noise\d+)?(_section\d+)?',
    well_stitching = '(' + '|'.join(wells) + ')(_subset\d+)?(_split\d+)?(_(cycle|)noise\d+)?',
    well_nonoise = '(' + '|'.join(wells) + ')(_subset\d+)?(_split\d+)?',
    well_nosplit = '(' + '|'.join(wells) + ')(_subset\d+)?',
    well_nosubset = '(' + '|'.join(wells) + ')',
    well_base = '(' + '|'.join(wells) + ')',

    tile = '\d\d\d\d',
    cycle = '|'.join(cycles_pt),

    path = '([^/]*/)*[^/.]*',
    path_nogrid = '((?!_grid\d)[^.])*',

    segmentation_type = '(cells|nuclei)(|bases)(|(expanded\d+))',


if type(config['phenotyping_channels'][0]) != list:
    config['phenotyping_channels'] = [config['phenotyping_channels']]

#assert len(set(config['sequencing_channels'])) == len(config['sequencing_channels'), (
    #"all channels in config['sequencing_channels' must be unique")

all_phenotyping_channels = []
for channels in config['phenotyping_channels']:
    #assert len(set(channels)) == len(channels) and all(chan not in all_phenotyping_channels for chan in channels), (
        #"all channels in config['phenotyping_channels'] must be unique")
    all_phenotyping_channels.extend(channels)

# Regex for use with wildcard constraints
phenotyping_channel_regex = '(' + '|'.join('{}|{}'.format(i, re.escape(name))
            for i, name in enumerate(all_phenotyping_channels)) + ')'
sequencing_channel_regex = '(' + '|'.join('{}|{}'.format(i, re.escape(name))
            for i, name in enumerate(config['sequencing_channels'])) + ')'
# only matches channels in both sequencing and phenotyping
any_channel_regex = ('(' + '|'.join(map(str, range(min(len(config['sequencing_channels']), len(config['sequencing_channels'])))))
    + '|' + '|'.join(map(re.escape, set(config['sequencing_channels']) & set(all_phenotyping_channels))) + ')')

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
    return '(' + ''.join('(_{}[^_]*)?'.format(name) for name in params) + ')'


    if kind is None and cycle is not None:
        kind = 'phenotyping' if cycle in phenotype_cycles else 'sequencing'
    if type(channel) == str:
        return config[kind+'_channels'].index(channel)
    return channel

def channel_index_phenotyping(channel):
    cycle = 0
    while channel not in config['phenotyping_channels'][cycle]:
        cycle += 1
    return cycle, config['phenotyping_channels'][cycle].index(channel)

def channel_index(channel, kind=None, cycle=None):
    if type(channel) == str:
        if kind is None and cycle is not None:
            kind = 'phenotyping' if cycle in phenotype_cycles else 'sequencing'

        if kind == 'phenotyping':
            cycle = phenotype_cycles.index(cycle)
            return config['phenotyping_channels'][cycle].index(channel)
        elif kind == 'sequencing':
            return config['sequencing_channels'].index(channel)
    return channel

ashlar_params_nooverlap  = [name.replace('_', '') for name in config['stitching']['ashlar'].keys()]
ashlar_params = ashlar_params_nooverlap + ['overlap', 'input', 'ashlar']
#ashlar_params = ['flip-x', 'flip-y', 'transpose', 'interp', 'filter-sigma', 'input']
#print (params_regex('channel', 'subpix', 'solver', *ashlar_params))

wildcard_constraints:
    ashlar_params = params_regex(*ashlar_params)


def find_input_file(wildcards=None, well=None, cycle=None):
    if well is None:
        well = [value for key, value in wildcards.items() if key[:4] == 'well'][0]
    if cycle is None:
        cycle = wildcards.cycle
    alternate_path = input_dir + '{well}/cycle{cycle}/raw.tif'.format(well=well, cycle=cycle)
    return config['inputfiles'].get(well, {}).get(cycle, alternate_path)


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
    print ('  segmentation directory:', segmentation_dir)
    print ('  phenotyping directory:', phenotyping_dir)
    print ('  output directory:', output_dir)
    print ()
    print ('Input found:')
    for well, files in inputfiles.items():
        print ('  ' + well + ':')
        print ('\n'.join('    {}: {}'.format(*pair) for pair in files.items()))
    print ('  Cycles:', ', '.join(cycles))
    print ()
    print ('  Wells:', ', '.join(wells))
    print ()

rule print_info:
    run:
        print_info()

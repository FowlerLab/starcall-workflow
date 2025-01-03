import os
import sys
import glob
import time

rawinput_dir = config.get('rawinput_dir', 'rawinput/')
log_dir = config.get('log_dir', 'logs/')

input_dir = config.get('input_dir', 'input/')
stitching_dir = config.get('stitching_dir', 'stitching/')
sequencing_dir = config.get('sequencing_dir', 'sequencing/')
phenotyping_dir = config.get('phenotyping_dir', 'phenotyping/')
output_dir = config.get('output_dir', 'output/')
qc_dir = config.get('qc_dir', output_dir + 'qc/')

stitching_input_dir = stitching_dir + input_dir
stitching_output_dir = stitching_dir + output_dir
sequencing_input_dir = sequencing_dir + input_dir
sequencing_output_dir = sequencing_dir + output_dir
phenotyping_input_dir = phenotyping_dir + input_dir
phenotyping_output_dir = phenotyping_dir + output_dir

phenotype_date = config.get('phenotype_date', 'phenotype')
#phenotype_cycle = config.get('phenotype_cycle', 'PT')
phenotype_cycles = config.get('phenotype_cycles', ['PT', 'PT1', 'PT2', 'PT3', 'PT4'])
phenotype_scale = config.get('phenotype_scale', 2)

dates = [] if not os.path.exists(rawinput_dir) else sorted(os.listdir(rawinput_dir))
dates_pt = dates.copy()
phenotype_dates = [date for date in dates if date[:len(phenotype_date)] == phenotype_date]
phenotype_cycles = phenotype_cycles[:len(phenotype_dates)]

if 'phenotype_channels' not in config:
    if phenotype_date in dates:
        import nd2
        phenotype_file = nd2.ND2File(glob.glob(rawinput_dir + '/phenotype/*.nd2')[0])
        phenotype_channels = phenotype_file.shape[1]
        phenotype_file.close()
        del phenotype_file
    else:
        phenotype_channels = 0
else:
    phenotype_channels = config['phenotype_channels']

if 'wells' not in config:
    #wells = sorted([path.replace('Well', '').partition('_')[0] for path in os.listdir(rawinput_dir + '/' + dates[0])])
    wells = sorted(list(set([path.partition('Well')[2].partition('_')[0] for path in glob.glob(rawinput_dir + '/*/*.nd2')])))
else:
    wells = config['wells']

for date in phenotype_dates:
    if date in dates_pt:
        dates.remove(date)

if 'cycles' not in config:
    cycles = ['{:02}'.format(i) for i in range(len(dates))]
else:
    cycles = config['cycles']

if phenotype_date in dates_pt:
    cycles_pt = cycles + phenotype_cycles
else:
    cycles_pt = cycles
#cycles_pt = sorted(cycles_pt)

cellpose_cyto_index = config.get('cellpose_cyto_index', 1)
cellpose_diameter = config.get('cellpose_diameter', 50)
cellpose_cycle = config.get('cellpose_cycle', cycles[-1] if len(cycles) else None)

alignment_channel = config.get('alignment_channel', 0)

apply_background_correction = config.get('apply_background_correction', False)

if 'subset' in config:
    tiles = tiles[::len(tiles)//10]

wildcard_constraints:
    well = '|'.join(wells),
    tile = '\d\d\d\d',
    cycle = '|'.join(cycles_pt),
    any_input_dir = '|'.join([input_dir, stitching_input_dir, sequencing_input_dir, phenotyping_input_dir]),
    any_output_dir = '|'.join([output_dir, stitching_output_dir, sequencing_output_dir, phenotyping_output_dir]),
    any_processing_dir = '|'.join([stitching_dir, sequencing_dir, phenotyping_dir]),
    filetype = '(\.[^/]+)|(/cycle(' + '|'.join(cycles_pt) + '))|',
    prefix = '(?!' + input_dir + ')(?!' + output_dir + ')([^/]*/)*[^/.]*',

def debug(*args, **kwargs):
    print (time.asctime() + ':', *args, **kwargs, file=sys.stderr)

def print_mem(name, mem_limit):
    import psutil
    import fisseq.utils
    #{key: fisseq.utils.human_readable(val) for key,val in psutil.Process().memory_info()._asdict().items()}
    mem = psutil.Process().memory_info().rss / 1000000
    debug ("Rule", name, "Limit", mem_limit, fisseq.utils.human_readable(mem_limit * 1000000),
                    "Using", mem, fisseq.utils.human_readable(mem * 1000000))

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

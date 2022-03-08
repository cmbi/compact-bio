# standard library imports
import argparse
import os
import sys
from importlib.util import spec_from_file_location,module_from_spec

# local imports
from rbomcl.main import main
from rbomcl import process_data as prd

def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="compact",
        description="CompaCt command line tool",
        )
    
    # arguments
    parser.add_argument('-v','--version', action='version',version=f'compact v0.0.1')
    parser.add_argument("settings",help="path of input settings file",
                        type=argparse.FileType('r'))
    parser.add_argument(
        "-p", type=float,default=0.9,
        help="rank biased overlap p parameter, default=0.9")
    parser.add_argument(
        "-m","--min-search-weight",type=float,default=0.999,
        help='minimum search weight for rank biased overlap, default=0.999'
        )
    parser.add_argument(
        '--th-criterium',choices=['percent','best'],default="percent",
        help="threshold criterium for determination of reciprocal top hit"
    )
    parser.add_argument(
        "--th-percent",type=int,default=1,
        help="if th-criterium='percent': top percentage to consider as top hits. default=1"
    )
    parser.add_argument(
        '--include_within',action="store_true",
        help='include within-sample interaction scores in clustering input'
    )
    parser.add_argument(
        '--wbratio',type=int,default=1,
        help="ratio of within/between score averages. default=1"
    )
    parser.add_argument(
        '-i','--mcl-inflation',default=2,
        help='inflation param of MCL tool, determines cluster granularity. default=2'
    )
    parser.add_argument(
        '-o','--output-loc',default=os.getcwd(),
        help='preferred location of the result folder. defaults to current working directory'
    )
    parser.add_argument(
        '-j','--job-name', default=None,
        help='name of current job, used in name of result folder name. defaults to collection names'
    )
    parser.add_argument(
        '--save-rthits',action="store_true",
        help="save reciprocal top hits to disk, in result folder"
    )
    parser.add_argument(
        '--report-threshold',type=float,default=0.5,
        help='fraction of samples threshold for reporting cluster membership'
    )
    parser.add_argument(
        '--perf-cluster-annot',action="store_true",
        help="perform automatic annotation of clusters using reference groups"
    )
    parser.add_argument(
        '--ref-file',type=argparse.FileType('r'),
        help='gmt format file with reference groups, when perf-cluster-annot is used'
    )
    parser.add_argument(
        '--ref-tag',default = None,
        help='tag of collection to use for automatic annotation.ids should match reference'
    )
    parser.add_argument(
        '--annot-ref-th',type=float,default=0.5,
        help='min fraction of reference group that must be present in cluster'
    )
    parser.add_argument(
        '--annot-mem-th',type=float,default=0.1,
        help='min fraction-present of clust members to be counted for annot-ref-th'
    )
    parser.add_argument(
        '-t','--processes',type=int,default=1,
        help="number of threads or processes to use. default=1"
    )

    args = parser.parse_args()

    ## some further processing/checking of parsed arguments##
    # close the settings input file wrapper
    args.settings.close()
    args.settings = args.settings.name
    if args.ref_file:
        args.ref_file.close()
        args.ref_file = args.ref_file.name
    # check result folder exists
    if not os.path.isdir(args.output_loc):
        raise ValueError(f'output location does not exist: {args.output_loc}')
    
    return args

def parse_settings(settings_fn):
    settings_name = os.path.basename(settings_fn).split('.')[0]
    spec = spec_from_file_location(settings_name,settings_fn)
    settings = module_from_spec(spec)
    spec.loader.exec_module(settings)
    
    sample_data = settings.sample_data

    mapping_data = settings.mapping_data

    return sample_data,mapping_data

def parse_profiles(fn_dict):
    """
    parse profiles for all complexomes,samples in fn_dict
    """
    sample_dict = {}
    for complexome,samplefiles in fn_dict.items():
        complexome_samples = {}
        for samplename,fname in samplefiles.items():
            sample = prd.parse_profile(fname)
            complexome_samples[samplename] = sample
        sample_dict[complexome] = complexome_samples
   
    return sample_dict

def parse_mappings(fn_dict):
    return {name:prd.parse_mapping(fn) 
            for name,fn in fn_dict.items()}

def get_nested_tags(corr_dict):
    """
    """
    nested_tags = {}
    for comp_name,comp_dict in corr_dict.items():
        nested_tags[comp_name] = list(comp_dict.keys())
    return nested_tags

def get_int_matrices(corr_dict):
    """
    """
    int_matrices = {}
    for comp_dict in corr_dict.values():
        int_matrices = {**int_matrices,**comp_dict}
    return int_matrices

def run():
    # parse arguments
    try:
        args = parse_arguments()
    except Exception as e:
        print(f'problem parsing arguments: {e}',file=sys.stderr)
        sys.exit()
    
    # parse input settings
    sample_data,mapping_data = parse_settings(args.settings)

    # parse data files
    try:
        samples = parse_profiles(sample_data)
    except Exception as e:
        print(f'problem parsing sample data: {e}',file=sys.stderr)
        sys.exit()

    try:
        mappings = parse_mappings(mapping_data)
    except Exception as e:
        print(f'problem parsing mapping data: {e}',file=sys.stderr)
        sys.exit()

    # run analysis
    try:
        #run analysis
        pass
    except Exception as e:
        print(f'problem during analysis: {e}',file=sys.stderr)

    return 0

if __name__ == "__main__":

    run()

    # print(args)

    # print("hello")
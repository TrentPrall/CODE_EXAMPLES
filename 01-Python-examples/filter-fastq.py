## modules ##
import os
from glob import glob
import sys

## parameters ##
# path to folder of ONT run folders, add run folders at sequencing continue to this directory
ont_run_folder_path = '/Volumes/TRENT_DATA3/01_ONT/Robot_Manuscript_Revision_3/SRA_Upload'

# docker invokation - save parameters that are used to invoke docker for this workflow
docker_run_cmd = 'docker run --user $(id -u):$(id -g) -it --read-only -v /slipstream:/slipstream -v $(pwd):/scratch -w /scratch dockerreg.chtc.wisc.edu/dhoconno/bulkvis:22510'

# ignore hidden files
def listdir_nohidden(path):
    out = []
    for f in os.listdir(path):
        if not f.startswith('.'):
            out.append(f)
        else:
            pass
    return out
            

def identify_available_folders(animal_ont_dir):
    '''for a given set of ONT folder data, create a data structure that holds the information necessary to run bulkvis.
    the data structure uses the ONT folder run name as the key and a tuple of the fastq_pass folder location and seq statistics as the value
    this allows the key to be used as the expansion unit in snakemake - each key gets expanded as a separate set of tasks that need to be run
    '''

    paths = glob(animal_ont_dir + '/*') # finds oc lab names
    AVAILABLE_RUNS = {}

    for i in paths:
        #print(i)
        # get folder name - use as dictionary key
        RUN_NAME = os.path.split(i)[-1]
        FLOW_CELL = listdir_nohidden(i)[0]
        #print('Flow Cell: ' + str(FLOW_CELL))

        # get path to fastq_pass files
        # there is an arbitrary folder between the run name and the fastq_pass
        FASTQ_PASS = glob(i + '/*/fastq_pass')
        # get path to adaptive sampling unblocked_read_ids.txt file
        if os.path.exists(i + '/' + FLOW_CELL + '/unblocked_read_ids.txt') == True:
            UNBLOCKED_READS = glob(i + '/*/unblocked_read_ids.txt')
        elif os.path.exists(i + '/' + FLOW_CELL + '/unblocked_read_ids.txt') == False:
            UNBLOCKED_READS = "NA"
        
        # make dictionary using RUN_NAME as key and FASTQ_PASS and UNBLOCKED_READS as list
        # only process runs where there are both FASTQ_PASS reads and UNBLOCKED_READS available
        if FASTQ_PASS and UNBLOCKED_READS:
            AVAILABLE_RUNS.update({RUN_NAME : [FASTQ_PASS[0], UNBLOCKED_READS[0]]})
    #print(AVAILABLE_RUNS)
    return(AVAILABLE_RUNS)

def check_file_permissions(animal_ont_dir):
    '''test permissions of files to process and make sure they are readable'''

    # loop through all files in folder tree
    # technically only need fastq files and sequencing stats
    # but probably a good idea to fix permissions if any files are unreadable

    unreadable_files = []

    for dirpath, dirs, files in os.walk(animal_ont_dir):
        for filename in files:
            fname = os.path.join(dirpath,filename)

            # collect files that aren't readable into a list
            if os.access(fname, os.R_OK) == False: unreadable_files.append(fname)
            
    # if there are any files that are unreadable, print an error message
    # list the unreadable files
    # then exit the script

    if unreadable_files:
        print('The following files cannot be read. This is almost due to permissions being set incorrectly. Change permissions with `chmod` and try again. ONLY CHANGE DIRECTORIES THAT CANNOT BE READ. IF YOU CHANGE THE ENTIRE PATH IT SNAKEMAKE WILL RERUN ALL FILES')
        print('\n'.join(unreadable_files))
        sys.exit('[insert sad trombone sound here]')

# specify folders to process
folders_to_process = identify_available_folders(ont_run_folder_path)
folder_names = list(folders_to_process.keys())
unblocked = {}
for keys in folders_to_process:
    unblocked[keys] = folders_to_process[keys][1]

# check that all files in folders_to_process are readable
check_file_permissions(ont_run_folder_path)

# docker cmd
docker_run_cmd = 'docker run --user $(id -u):$(id -g) -v $(pwd):/scratch -w /scratch pegi3s/seqkit grep '

'''
TEST
docker run --user $(id -u):$(id -g) -v $(pwd):/scratch -w /scratch pegi3s/seqkit grep \
--pattern-file 20210212_25181_r02072_N4/unblocked_reads.txt \
20210212_25181_r02072_N4/20210212_25181_r02072_N4-merged_reads.fq.gz
'''
